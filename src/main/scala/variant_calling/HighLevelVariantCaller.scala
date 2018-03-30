package variant_calling

/**
  * Author: Alex N. Salazar
  * Created on 26-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

import java.io.{File, PrintWriter}

import utilities.GFAutils.ConstructGFA
import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile}

import scala.collection.immutable.HashSet
import quiver._
import utilities.HLGG

import scala.annotation.tailrec

object HighLevelVariantCaller extends ConstructGFA {

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
      highLevelVariantCaller(config)
    }
  }

  def highLevelVariantCaller(config: Config): Unit = {
    //load HLGG from given GFA file
    val hlgg = new HLGG(config.hlgg)
    println("Using " + config.reference + " as reference genome architecture" + timeStamp)
    //extract map containing reference sequence identifier and the corresponding arrows in positional order
    val ref_arrows = hlgg.sequence_set.filter(_._1.contains(config.reference))
    assert(!ref_arrows.isEmpty, "Could not find reference genome containing the follownig substring: " + config
      .reference)
    //create a hashset containig all nodes element of the reference genome
    val ref_nodes = ref_arrows.flatMap(_._2).foldLeft(HashSet[Int]())((set, node) => set + node)
    //create hashset containing all sequence identifiers of the reference
    val ref_labels = ref_arrows.map(_._1).foldLeft(HashSet[String]())((set, label) => set + label)
    println("Reference has " + ref_arrows.map(_._2.size).sum + " arrows with " + ref_nodes.size + " nodes and " +
      ref_labels.size + " individual labels" + timeStamp)


    /**
      * Function to write quiver graph into GFA format
      *
      * @return
      */
    def graphToGFA: (hlgg.HLGG, PrintWriter) => Unit = (graph, pw) => {
      graph.labNodes.foreach(node => pw.println(constructSegmentLine(node.vertex)))
      graph.labEdges.foreach(edge => pw.println(constructLinkLine(edge.from, edge.to)))
    }

    /**
      * Function to convert a path (e.g. maximally-labelled path) encoded as a GFA line
      *
      * @return
      */
    def pathToGFA: (Int, List[LEdge[Int, Set[String]]], PrintWriter) => Unit = (id, path, pw) => {
      val t = getPath(path)
      val label = path.map(_.label).distinct
      assert(label.size == 1)
      pw.println(Seq("P", id, t.mkString("+,"), label.head.mkString(",")).mkString("\t"))
    }


    /**
      * Function to check whether a given label contains a specified reference genome
      *
      * @return
      */
    def hasRef: Set[String] => Boolean = label => label.exists(ref_labels(_))

    /**
      * Function to obtain the corresponding path of nodes given a corresponding path of arrows.
      *
      * @return
      */
    def pathComposition(subgraph: hlgg.HLGG): Seq[LEdge[Int, Set[String]]] => Seq[Int] = path => {
      //get size of path
      val size = path.size
      //iterate through path joined with index
      path.zipWithIndex.foldLeft(Seq[Int]()) { case (composition, (arrow, index)) =>
        //add only the tail node of an arrow if it's not the last one
        if (index + 1 != path.size) composition.:+(arrow.from)
        else {
          //add only the tail node if the head node is not a leaf in the current subgraph
          if (!subgraph.isLeaf(arrow.to)) composition.:+(arrow.from) else (composition.:+(arrow.from)).:+(arrow.to)
        }
      }
    }

    def getPath: List[LEdge[Int, Set[String]]] => Seq[Int] = path => {
      val size = path.size
      path.zipWithIndex.foldLeft(Seq[Int]()) { case (composition, edge) => {
        if (edge._2 + 1 == size) composition.:+(edge._1.from).:+(edge._1.to) else composition.:+(edge._1.from)
      }
      }
    }

    /**
      * Method to identify and return all connected components in a quiver graph. Starting from some root, performs
      * BFS-traversal to identify all nodes in the connected component. Then, removes them and returns a new subgraph
      * . Methods continues until graph is empty.
      *
      * @param remaining_subgraph
      * @param subgraph
      * @return
      */
    @tailrec def getCCs(remaining_subgraph: hlgg.HLGG,
                        subgraph: List[hlgg.HLGG],
                        collectSingletons: Boolean): List[hlgg.HLGG] = {
      /**
        * Tail-recursive method to traverse through the reference sequence path given it's root node and
        * obtain all nodes in the path
        *
        * @param nodesToVisit
        * @param visited
        * @return
        */
      @tailrec def undirectedBfs(nodesToVisit: Seq[Int],
                                 visited: HashSet[Int]): HashSet[Int] = {
        if (nodesToVisit.isEmpty) visited
        else {
          //get current node
          val current_node = nodesToVisit.head
          //get the context of the current node
          val current_context = remaining_subgraph.context(current_node)
          //get all valid neighbours
          val neighbours = (current_context.outEdges ++ current_context.inEdges).flatMap(x => List(x.from, x.to))
            .filterNot(x => x == current_node || visited(x))
          //continue
          undirectedBfs(nodesToVisit.tail ++ neighbours, (visited + current_node) ++ neighbours)
        }
      }

      if (remaining_subgraph.isEmpty) subgraph
      else {
        //find some random root
        val root = remaining_subgraph.roots.head
        //perform undirected bfs-traversal
        val cc = undirectedBfs(Seq(root), HashSet(root))
        //remove nodes from global graph in the connected component
        val updated_graph = cc.foldLeft(remaining_subgraph)((b, a) => b.removeNode(a))
        //get subgraph of only the conncted component
        val cc_subgraph = hlgg.empty_HLGG.addNodes(remaining_subgraph.labNodes.filter(x => cc(x.vertex)))
          .addEdges(remaining_subgraph.labEdges.filter(x => cc(x.from) && cc(x.to)))
        //act accordingly based on size of connected component
        if (collectSingletons) getCCs(updated_graph, cc_subgraph :: subgraph, collectSingletons)
        else {
          cc.size match {
            case 1 => getCCs(updated_graph, subgraph, collectSingletons)
            case _ => getCCs(updated_graph, cc_subgraph :: subgraph, collectSingletons)
          }
        }
      }
    }

    /**
      * Tail-recursive method to obtain the maximum-labeled path in a quiver graph. Starting from some root and it's
      * label, the maximum-labeled path is the longest traversal of arrows with identical labels.
      *
      * @param subgraph
      * @param maximal_paths
      * @return
      */
    @tailrec def maximalLabeledPath(subgraph: hlgg.HLGG,
                                    maximal_paths: List[List[LEdge[Int, Set[String]]]]): List[List[LEdge[Int, Set[String]]]] = {

      /**
        * Function to obtain the maximal-labeled path of from the current subgraph given a starting node and a
        * label using DFS.
        *
        * @return
        */
      def traverse: (Int, Set[String]) => Vector[LEdge[Int, Set[String]]] = (root, label) => {
        subgraph.xdfsWith(Seq(root),
          (c: hlgg.HLGC) => c.outEdges.filter(_.label == label).map(_.to),
          (c: hlgg.HLGC) => c.outEdges.filter(_.label == label)
        ).flatten
      }

      //get roots of the graph
      val roots = {
        val r = subgraph.roots.toList
        if (!r.isEmpty) r else subgraph.nodes.take(1)
      }
      //if there are no roots, we have finished traversing the graph
      if (roots.isEmpty) maximal_paths
      else {
        //find singletons, disconnected components made of one node. by definition, these are roots and leafs nodes
        val singletons = roots.filter(subgraph.isLeaf(_))
        //remove all singletons first, if they exist
        if (!singletons.isEmpty) maximalLabeledPath(subgraph.removeNodes(singletons), maximal_paths)
        else {
          //get all outgoing root arrows; note: there can be multiple roots
          val root_edges = roots.flatMap(x => subgraph.context(x).outEdges).toList
          //for each outgoint arrow, perform maximum labeled traversal based on DFS
          val traversals = root_edges.map(x => traverse(x.from, x.label).toList)
          //remove all arrows in the traversals from the current graph and add traversals to maximal_paths collection
          maximalLabeledPath(traversals.flatten.foldLeft(subgraph)((s, e) => s.removeLEdge(e)),
            traversals.foldLeft(maximal_paths)((b, a) => a :: b))
        }
      }
    }

    println("Getting connected components" + timeStamp)
    //get all resulting disconnect components
    val ccs = getCCs(hlgg.graph, List(), true)
    println("--Found " + ccs.size + " connected components" + timeStamp)
    println("Performing maximally-labelled path traversals:" + timeStamp)
    //iterate through each disconnected component and identify structural variant calls
    ccs.foldLeft(0)((cc_id, cc) => {
      println("--Current connected component " + cc_id + " has " + cc.labEdges.size + " total edges" + timeStamp)
      //remove all arrows in the reference direction; as well as singleton cycles (e.g. annotation errors)
      val edges_to_remove = cc.labEdges.filter(x => hasRef(x.label) || x.from == x.to)
      println("--Removing " + edges_to_remove.size + " reference edges" + timeStamp)
      val ref_edges_removed = edges_to_remove.foldLeft(cc)((b, a) => b.removeLEdge(a))
      //get family of subgraphs
      val family_subgraphs = getCCs(ref_edges_removed, List(), false)
      println("--Partitioned into a family of subgraphs of " + family_subgraphs.size + " subgraphs" + timeStamp)
      //create output file
      val output = new PrintWriter(config.outputDir + "/connected_component." + cc_id + ".gfa")
      //output family of subgraphs to gfa
      family_subgraphs.foreach(s => graphToGFA(s, output))
      println("--Computing maximially-labelled paths" + timeStamp)
      //iterate throug each subgraph and characterize maximally-labelled paths
      family_subgraphs.foldLeft(0) {
        case (local_id, family_subgraph) => {
          //get mlp for current family subgraph
          val maximially_labelled_paths = maximalLabeledPath(family_subgraph, List())
          //if(family_subgraph.contains(70953)) maximially_labelled_paths.foreach(println)
          //iterate through each mlp and outpu to gfa
          maximially_labelled_paths.zipWithIndex.foreach {
            case (mlp, index) => pathToGFA(index + 1 + local_id, mlp, output)
          }
          (local_id + maximially_labelled_paths.size)
        }
      }
      output.close
      cc_id + 1
    })
  }
}
