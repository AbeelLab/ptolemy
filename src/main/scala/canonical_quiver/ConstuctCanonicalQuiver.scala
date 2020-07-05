package canonical_quiver

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory}
import utilities.GFAutils.{ConstructGFA}
import utilities.ConfigHandling

import scala.collection.immutable.{HashMap, HashSet}

/**
  * Author: Alex N. Salazar
  * Created on 16-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ConstuctCanonicalQuiver extends ConstructGFA {

//  case class Config(
//                     syntenicAnchors: File = null,
//                     database: File = null,
//                     outputDir: File = null,
//                     verbose: Boolean = false,
//                     msa: Boolean = false,
//                     isCircular: Boolean = false,
//                     dump: Boolean = false
//                   )

  def main(args: Array[String]) {
    val defaultValues = ConfigHandling.fullConfig()
    val parser = new scopt.OptionParser[ConfigHandling.fullConfig]("canonical-quiver") {
      opt[File]('s', "syntenic-anchors") required() action { (x, c) =>
        c.copy(syntenicAnchors = x)
      } text ("Path to file containing syntenic anchors")
      opt[File]('d',"db") required() action { (x, c) =>
        c.copy(database = x)
      } text ("\n"+" "*27+"Directory path of database (i.e. output-directory used "+
      "\n"+" "*27+"for the 'extract' module)")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory")
      note("\nOPTIONAL FLAGS")
//      opt[Unit]("circular") action { (x, c) =>
//        c.copy(isCircular = true)
//      } text ("Genome can be circular (default is false)")
      opt[Unit]("msa-groups") action { (x, c) =>
        c.copy(msa = true)
      } text ("\n"+" "*27+"Turn on to output file of syntenic anchors to database"+
        "\n"+" "*27+"for inducing Multiple-Sequence-Alignment in each"+
        "\n"+" "*27+"syntenic anchor ("+defaultValues.msa+" by default)")
      opt[Unit]("dump") action { (x, c) =>
        c.copy(dump = true)
      } text ("\n"+" "*27+"Write intermediate tables to disk ("+defaultValues.dump+" by default)")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      } text("\n"+" "*27+"Display extra process information (default is "+ defaultValues.verbose+")")
    }
    parser.parse(args, ConfigHandling.fullConfig()).map { parsedConfig =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(parsedConfig.database)
      // handle options and flags for the current module
      val config = ConfigHandling.parameterManager(parsedConfig, "canonical-quiver")
      constructHLGG(config)
    }
  }

  def constructHLGG(config: ConfigHandling.fullConfig): Unit = {
    println(timeStamp + "Fetching hashmap Z")
    val path_hashmap_Z = config.database.listFiles().find(_.getName == "global_z.txt").get
    val path_hashmap_Y = config.database.listFiles().find(_.getName == "global_y.txt").get
    val path_orfids = config.database.listFiles().find(_.getName == "orf2id_mapping.txt").get
    println(timeStamp + "Fetching starting node ID")
    //get last orf id to be used as the starting node id for new nodes
    val node_id = {
      //open file as iterator
      val iterator = openFileWithIterator(path_orfids)
      //get last orf id, increment +1
      iterator.foldLeft(0)((last, line) => {
        val id = line.split("\t").last.toInt
        if (id > last) id else last
      }) + 1
    }
    //get total number of genomes
    val total_genomes = openFileWithIterator(path_hashmap_Y).toList.size
    println(timeStamp + "--Starting node is: " + node_id)
    println(timeStamp + "Formatting syntenic anchors")
    //open syntenic anchors and format to hashmap as original orf id -> (assigned node id, genome count)
    val (syntenic_anchors, node_to_genome_count, next_id) =
      openFileWithIterator(config.syntenicAnchors).foldLeft(HashSet[Set[Int]]())((sa, line) => {
      val anchors = {
        val tmp = line.split("\t")
        (tmp.take(1) ++ tmp(1).split(",")).map(_.toInt).toSet
      }
      sa + (anchors)
    }).foldLeft((HashMap[Int,Int](), HashMap[Int,Int](), node_id)){ case((map, counts, id), sa) => {
        if(sa.size > total_genomes){
          println("--WARNING: unexpected number of nodes in a syntenic anchor. Discarding: " + sa.mkString(","))
          (map, counts, id)
        } else {
          (sa.foldLeft((map))((local_map, orf) => local_map + (orf -> id)), counts + (id -> sa.size), id + 1)
        }
      }}

    /**
      * Function to fetch node ID given an ORF ID
      * @return node ID as INT
      */
    def getNodeID: Int => Int = orf => {
      val get = syntenic_anchors.get(orf)
      if (get == None) orf else get.get
    }

    if (config.dump) {
      val pw = new PrintWriter(config.outputDir + "/node2genome_count.txt")
      node_to_genome_count.foreach(x => pw.println(x._1 + "\t" + x._2))
      pw.close
    }

    if(config.msa){
      println(timeStamp + "User specified MSA groups to database. Writing to disk")
      val pw = new PrintWriter(config.database + "/msa_groups.txt")
      syntenic_anchors.toList.groupBy(_._2).foreach(sa => pw.println(sa._1 + "\t" + sa._2.map(_._1).mkString(",")))
      pw.close
    }

    println(timeStamp + "Constructing canonical quiver")
    //iterate through each sequence and construct HLGG in context of syntenic anchors
    val (hlgg_edges, hlgg_nodes) =
      openFileWithIterator(path_hashmap_Z)
        //iterate through each sequence with global adjacency map and node set
        .foldLeft((HashMap[Int, Set[Int]](), HashSet[Int]())) { case ((adj_map, node_set), _sequence) => {
        //get sequence id and sequence of orf ids
        val (sequence, orfs) = getSequence(_sequence)
        orfs.size match {
          case 1 => (adj_map, node_set + orfs.head)
          case _ => {
            //create the sliding window depending on whether genomes are cicular or not
            val orf_slide = {
              val tmp = orfs.sliding(2)
              if(!config.isCircular) tmp else tmp ++ Iterator(Seq(orfs.last, orfs.head))
            }
            //iterate through each orf as (node, edge)
            orf_slide.foldLeft((adj_map, node_set)) { case ((local_adj_map, local_node_set), adj) => {
              //determine node and edge id
              val node_id = getNodeID(adj(0))
              val edge_id = getNodeID(adj(1))
              //get current edges for current node
              val current = local_adj_map.getOrElse(node_id, Set[Int]())
              //update accordingly
              (local_adj_map + (node_id -> current.+(edge_id)), (local_node_set + (node_id)) + edge_id)
            }
            }
          }
        }
      }
      }

    println(timeStamp + "Constructed canonical quiver with " + hlgg_nodes.size + " nodes and " + hlgg_edges.map(_._2.size).sum +
      " edges")
    val path2Nodes = scala.collection.mutable.Map[String,Seq[Boolean]]()
    println(timeStamp + "Writing GFA file to disk")
    val sortedNodeID = hlgg_nodes.toSeq.sorted
    val pw = new PrintWriter(config.outputDir + "/canonical_quiver.gfa")
    pw.println(getGFAHeader)
    hlgg_nodes.foreach(x => pw.println(constructSegmentLine(x) + addGenomeCountField(node_to_genome_count.get(x))))
    hlgg_edges.foreach { case (node, edges) => edges.foreach(edge => pw.println(constructLinkLine(node, edge))) }
    //create paths for each sequence
    openFileWithIterator(path_hashmap_Z).foreach(entry => {
      //parse
      val (sequence, orfs) = getSequence(entry)
      //exchange orf id to new assigned node id if it was involved in an anchor
      val orf_to_node_ids = orfs.map(orf => {
        val get = syntenic_anchors.get(orf)
        if (get == None) orf + "+" else get.get + "+"
      })
      //output paths
      pw.println(Seq("P", sequence, orf_to_node_ids.mkString(","), "*").mkString("\t"))

      // boolean (falses) list for the matrix table for current path
      val nodeInPath = scala.collection.mutable.Seq.fill(hlgg_nodes.size)(false)
      // loop over the nodes of the path
      for (node <- orf_to_node_ids) {
        // change the current node label of sorted sequence of node ID to true
        nodeInPath.update(sortedNodeID.indexOf(node.dropRight(1).toInt), true)
      }
      // assign the current nodeInPath sequence to the map
      path2Nodes(sequence) = nodeInPath
    })

    // map to store the relation between strains and paths
    val strain2Paths = scala.collection.mutable.Map[String, List[String]]()
    // save the strains into the file
    openFileWithIterator(path_hashmap_Y).foreach(line => {
      val split = line.split("\t")
      pw.println(Seq("G", split.head, split(1)).mkString("\t"))
      strain2Paths(split.head) = split(1).split(",").toList
    })
    pw.close

    val pw_orfids = new PrintWriter(config.database + "/orf2node_id.txt")
    openFileWithIterator(path_orfids).foreach(line => {
      val split = line.split("\t")
      val orf_id = split(2).toInt
      val node_id = getNodeID(orf_id)
      pw_orfids.println(Seq(split(0), split(1), orf_id, node_id).mkString("\t"))
    })
    pw_orfids.close

    // from strain2Paths and path2Nodes, strain2Nodes shoould be computed with 
    // a logical OR between the path nodes of each strain
    val strain2Nodes = scala.collection.mutable.Map[String, Seq[Boolean]]()
    scala.collection.immutable.ListMap(strain2Paths.toSeq.sortBy(_._1):_*).foreach{
      case (strain_i, pathList) =>
        var nodeList = scala.collection.mutable.Seq.fill(hlgg_nodes.size)(false)
        pathList.foreach{case path_i =>
          for (i <- nodeList.indices)
            nodeList(i) = nodeList(i) || path2Nodes(path_i)(i)
        }
        strain2Nodes(strain_i) = nodeList
    }

    // print the matrix table for the canonical quiver
    val pws = new PrintWriter(config.outputDir + "/canonical_quiver.strains.out")
    // print the header
    pws.print("Node\t")
    strain2Nodes.foreach{p => pws.print(p._1 + "\t")}
    pws.println("")
    for ((nodeID, i) <- sortedNodeID.zipWithIndex) {
      pws.print(nodeID + "\t")
      strain2Nodes.foreach{p => if (p._2(i)) pws.print("1\t") else pws.print("0\t")  }
      pws.println("")
    }
    pws.close

    println(timeStamp + "Successfully completed!")
  }

  /**
    * Function to parse file containing hashmap Z
    *
    * @return Tuple (String, Seq[Int]) representing Sequence ID and sorted ORFs
    */
  def getSequence: String => (String, Seq[Int]) = line => {
    val split = line.split("\t")
    (split.head, split(1).split(",").map(_.toInt).toSeq)
  }

}
