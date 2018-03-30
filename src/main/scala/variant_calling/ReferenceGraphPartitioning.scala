package variant_calling

import java.io.{File, PrintWriter}

import utilities.FileHandling.{verifyDirectory, verifyFile, timeStamp}
import utilities.GFAutils.ConstructGFA
import utilities.HLGG

import scala.annotation.tailrec
import scala.collection.immutable.HashSet

/**
  * Author: Alex N. Salazar
  * Created on 21-3-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */



object ReferenceGraphPartitioning extends ConstructGFA {

  case class Config(
                     hlgg: File = null,
                     outputDir: File = null,
                     reference: String = null,
                     verbose: Boolean = false,
                     dump: Boolean = false
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("construct-hlgg") {
      opt[File]('h', "hlgg-gfa") required() action { (x, c) =>
        c.copy(hlgg = x)
      } text ("GFA file of the HLGG.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory.")
      opt[String]("reference") required() action { (x, c) =>
        c.copy(reference = x)
      } text ("Reference genome. This is used to define a 'natural-flow' of a quiver which is then used to identify " +
        "all structural variants that deviate from it.")
      note("\nOPTIONAL\n")
      opt[Unit]("dump") action { (x, c) =>
        c.copy(dump = true)
      }
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      }

    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.hlgg)
      partitionGraph(config)
    }
  }

    def partitionGraph(config: Config): Unit = {
      println("Loading canonical quiver" + timeStamp)
      //load HLGG from given GFA file
      val hlgg = new HLGG(config.hlgg)
      //extract map containing reference sequence identifier and the corresponding arrows in positional order
      val ref_arrows = hlgg.sequence_set.filter(_._1.contains(config.reference))
      assert(!ref_arrows.isEmpty, "Could not find reference genome containing the follownig substring: " + config
        .reference)
      //create a hashset containig all nodes element of the reference genome
      val ref_nodes = ref_arrows.flatMap(_._2).foldLeft(HashSet[Int]())((set, node) => set + node)
      //create hashset containing all sequence identifiers of the reference
      val ref_labels = ref_arrows.map(_._1).foldLeft(HashSet[String]())((set, label) => set + label)
      println("Using " + config.reference + " genome architecture as reference" + timeStamp)

      /**
        * Method to identify and return all connected components in a quiver graph. Starting from some root, performs
        * BFS-traversal to identify all nodes in the connected component. Then, removes them and returns a new subgraph
        * . Methods continues until graph is empty.
        *
        * @param remaining_subgraph
        * @param subgraph
        * @return
        */
      @tailrec def getCCs(remaining_subgraph: hlgg.HLGG, subgraph: List[hlgg.HLGG]): List[hlgg.HLGG] = {
        if (remaining_subgraph.isEmpty) subgraph
        else {
          //find some random root
          val root = remaining_subgraph.roots.head
          //perform bfs-traversal
          val cc = remaining_subgraph.udfs(Seq(root)).toList
          //act accordingly based on size of connected component
          cc.size match {
            //if the cc is composed of a singleton node, remove it and move on
            case 1 => getCCs(remaining_subgraph.removeNodes(cc), subgraph)
            //else, extract subgraph from traversal and add as connected component
            case _ => getCCs(remaining_subgraph.removeNodes(cc), remaining_subgraph.subgraph(cc) :: subgraph)
          }
        }
      }
      println("Traversing canonical quiver and obtaining connected components" + timeStamp)
      //get all connected components in the canonical quiver
      val all_ccs = getCCs(hlgg.graph, List())
      println("Found " + all_ccs.size + " connected components" + timeStamp)

      /**
        * Function to traverse through each connected component (CC) and partition the CC based on the chose
        * reference sequence path
        * @param ccs
        * @param colours
        * @return
        */
      @tailrec def partitionCCs(ccs: List[hlgg.HLGG], colours: List[Set[Int]]): List[Set[Int]] = {
        ccs match {
            //return sets of reference/non-reference sequences
          case Nil => colours
          case cc :: tail => {
            //get all reference labels that exist in the current connected component, if any
            val current_reference_labels = cc.labNodes.foldLeft(Set[String]())((rlabels, node) => {
              val filtered_label = node.label.filter(ref_labels(_))
              if(filtered_label.isEmpty) rlabels else rlabels ++ filtered_label
            })
            //current cc contains sequences not present in the reference
            if(current_reference_labels.isEmpty) {
              println("--The current connected component does not have a reference sequence" + timeStamp)
              partitionCCs(tail, cc.nodes.toSet :: colours)
            }
            else {
              println("--The current connected component has " + current_reference_labels.size +
                " reference sequences" + timeStamp)
              //iterate through each reference label and extract all nodes part of the reference sequence traversal
              val colours_updated = current_reference_labels.foldLeft(colours){case (local_colours,ref) => {
                //get root of reference sequence
                val root = ref_arrows(ref).head
                //get all other reference labels to NOT traverse
                val anti_labels = current_reference_labels.foldLeft(HashSet[String]())((b,a) => if(a == ref) b else b + a)
                println("----Processing " + ref)
                /**
                  * Tail-recursive method to traverse through the reference sequence path given it's root node and
                  * obtain all nodes in the path
                  * @param nodesToVisit
                  * @param visited
                  * @return
                  */
                @tailrec def undirectedBfsWith(nodesToVisit: Seq[Int],
                                               refContained: Set[Int],
                                               visited: HashSet[Int]): (Set[Int], HashSet[Int]) = {
                  if(nodesToVisit.isEmpty) (refContained, visited)
                  else {
                    //get current node
                    val current_node = nodesToVisit.head
                    //get the context of the current node
                    val current_context = cc.context(current_node)
                    //get all valid neighbours
                    val neighbours =
                      (current_context.outEdges ++ current_context.inEdges).map(_.to).filterNot(x => x == current_node ||
                      visited(x) || cc.context(x).label.exists(anti_labels(_)))
                    //continue
                    undirectedBfsWith(nodesToVisit.tail ++ neighbours,
                      (if(current_context.label.exists(_ == ref)) refContained + current_node else refContained),
                      (visited + current_node) ++ neighbours)
                  }
                }
                val ref_nodes = cc.labNodes.foldLeft(HashSet[Int]())((b,a) =>
                  if(a.label.exists(_ == ref)) b + a.vertex else b)
                val (ref_nodes_only, all_visited) = undirectedBfsWith(ref_nodes.toSeq, Set(), ref_nodes)
                println("----Found " + ref_nodes_only.size + " reference nodes")
                (ref_nodes_only :: local_colours)
              }}
              partitionCCs(tail, colours_updated)
            }
          }
        }
      }
      println("Iterating through each connected component and performing reference-guide graph partitioning" + timeStamp)
      val partitioned_graph = partitionCCs(all_ccs, List())
      println("Writing partitions to disk" + timeStamp)
      val pw = new PrintWriter(config.outputDir + "/reference_partitioned.txt")
      partitioned_graph.foreach(x => pw.println(x.mkString(",")))
      pw.close
      println("Successfully completed!" + timeStamp)
    }

}
