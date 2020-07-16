package alignment

/**
  * Author: Alex N. Salazar
  * Created on 15-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */


import java.io.{File, PrintWriter}

import atk.FastaIterator
import utilities.FileHandling.{timeStamp, verifyDirectory, verifyFile, getFileName, getParentDirectory}
import utilities.GFAutils.ReadGFA
import utilities.DatabaseUtils.PtolemyDB
import utilities.MinimizerUtils.MMethods
import utilities.ConfigHandling

object BuildIndex extends ReadGFA with PtolemyDB with MMethods {

//  case class Config(
//                     canonicalQuiver: File = null,
//                     database: File = null,
//                     kmerSize: Int = 15,
//                     windowSize: Int = 3,
//                     // maxNodeFreq: Double = 0.01, // unused variable!
//                     verbose: Boolean = true,
//                   )

  def main(args: Array[String]) {
    val defaultValues = ConfigHandling.fullConfig()
    val parser = new scopt.OptionParser[ConfigHandling.fullConfig]("build-index") {
      opt[File]('c', "canonical-quiver") required() action { (x, c) =>
        c.copy(canonicalQuiver = x)
      } text ("Path to canonical quiver in GFA-format")
      opt[File]('d',"db") required() action { (x, c) =>
        c.copy(database = x)
      } text ("\n"+" "*27+"Directory path of database (i.e. output-directory used "+
      "\n"+" "*27+"for the 'extract' module)")
      note("\nOPTIONAL FLAGS")
      opt[Int]('k', "kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("\n"+" "*27+"Kmer sizer used during alignment (default is "+defaultValues.kmerSize+")")
      opt[Int]('w', "window-size") action { (x, c) =>
        c.copy(windowSize = x)
      } text ("Minimizer window size (default is "+defaultValues.minimizerWindow+")")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      } text("\n"+" "*27+"Display extra process information (default is "+ defaultValues.verbose+")")
    }
    parser.parse(args, ConfigHandling.fullConfig()).map { parsedConfig =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(parsedConfig.database)
      // handle options and flags for the current module
      val config = ConfigHandling.parameterManager(parsedConfig, "index-graph")
      //check whether canonical quiver exists prior to runing the indexing
      verifyFile(config.canonicalQuiver)
      buildIndex(config)
    }
  }

  def buildIndex(config: ConfigHandling.fullConfig): Unit = {
    println(timeStamp + "Loading canonical quiver")
    //load all nodes in the canonical quiver
    val gfa_nodes = loadNodesGFA(config.canonicalQuiver)
    println(timeStamp + "Found " + gfa_nodes.size + " nodes in the canonical quiver")
    println(timeStamp + "Loading ORF ID to node ID schema")
    //create map from orf ID to node IDs
    val orfid2nodeid = {
      val tmp = loadORFid2Nodeid(config.database)
      //get all found node IDs
      val found_nodeids = tmp.map(_._2).toSet
      //sanity checks
      gfa_nodes.foreach(x => assert(found_nodeids(x), "Could not find mapping for canonical quiver node: " + x))
      found_nodeids.foreach(x => assert(gfa_nodes(x), "Could not find mapping for database node: " + x))
      tmp
    }
    //call the curried method for fetching minimizers first setting standard parameters
    val generalMinimizerFetcher = fetchMinimizers(config.kmerSize, config.windowSize, config.verbose) _

    println(timeStamp + "Creating global kmer index for nodes")
    //create a global kmer index: hashtable with kmer as keys and a map of node id -> start position as values
    // as well as global node info index: hashtable with node id as keys and a map of
    // orf id -> (avg minimizer setsize, avg sequence size)
    val (global_kmer_index, global_nodeinfo) = {
      //fetch all ORF sequences files
      config.database.listFiles().filter(_.getName.endsWith("orfs.sequences.fasta"))
        //iterate through each one and update kmer index (Map[kmer -> Map[Node, starting positions]])
        .foldLeft((Map[Int, Map[Int, Set[(Int, Int)]]](), Map[Int, Map[Int, (Int, Int)]]())) {
        case ((kmer_index, nodeinfo), _fasta) =>
          //call for fetching minimizers from ORFs
          generalMinimizerFetcher(Option(orfid2nodeid), Option(nodeinfo))(_fasta, kmer_index)
      }
    }

    println(timeStamp + "Creating global node index")
    //using the global kmer index, create hashtable with node ids as keys and ordered sequence of starting positions
    val global_node_index = {
      //iterate through each kmer and it's map of node id's, ultimately updating the node index
      global_kmer_index.foldLeft(Map[Int, Set[(Int, Int, Int)]]()) { case (node_index, (kmer, node_map)) => {
        //iterate through each node in the node map and update the global starting positions
        node_map.foldLeft(node_index) { case (local_node_index, (node_id, starting_positions)) => {
          //get current starting positions
          val current_starting_positions = local_node_index.getOrElse(node_id, Set[(Int, Int, Int)]())
          //update starting positions
          local_node_index + (node_id -> (current_starting_positions ++ starting_positions.map(x => (kmer, x._1, x._2))))
        }}
      }}
        //sort starting positions in ascending order
        .mapValues(_.toSeq.sortBy(x => (x._2, x._1)))
    }

    println(timeStamp + "Creating global kmer index for intergenic sequences")
    //creartea global kmer index but for intergenic sequences
    val global_inter_kmer_index = {
      val tmp = config.database.listFiles().filter(_.getName.endsWith("intergenic.sequences.fasta"))
        //iterate through each one and update kmer index (Map[kmer -> Map[Node, starting positions]])
        .foldLeft((Map[Int, Map[Int, Set[(Int, Int)]]](), Map[Int, Map[Int, (Int, Int)]]())) {
        //call for fetching minimizers from intergenic sequences
        case ((kmer_index, nodeinfo), _fasta) => generalMinimizerFetcher(None, None)(_fasta, kmer_index)
      }
      tmp._1
    }

    println(timeStamp + "Writing to disk")
    //fetch name of canonical quiver file
    val name = getFileName(config.canonicalQuiver)
    //fetch parent directory
    val parent_directory = getParentDirectory(config.canonicalQuiver)
    //create output files for database indeces
    val pw_nki = new PrintWriter(parent_directory + "/" + name + ".nki")
    val pw_nki_distribution = new PrintWriter(parent_directory + "/" + name + ".nkid")
    val pw_iki = new PrintWriter(parent_directory + "/" + name + ".iki")
    val pw_ni = new PrintWriter(parent_directory + "/" + name + ".ni")
    val pw_nii = new PrintWriter(parent_directory + "/" + name + ".nii")
    //output kmer index
    // this is different the data structure constructed above, it's structured: <kmer hash>\t<node id>,<node id>...
    global_kmer_index.foreach { case (kmer, node_map) => {
      val local_nodes = node_map.foldLeft(Set[Int]())((b, a) => b + (a._1))
      pw_nki.println(kmer + "\t" + local_nodes.mkString(","))
      pw_nki_distribution.println(kmer + "\t" + local_nodes.size)
    }}
    //output intergenic kmer index
    global_inter_kmer_index.foreach{case (kmer, empty_map) => pw_iki.println(kmer)}
    //output node index
    global_node_index.foreach { case (node, starting_positions) =>
      pw_ni.println(node + ":" + starting_positions.map(x => (x._1 + "_" + x._2 + "_" + x._3)).mkString(","))
    }
    //summarize node info as: node id -> (average minimizer set size, average sequence size)
    global_nodeinfo.mapValues(orf_ids => (mean(orf_ids.map(_._2._1)), mean(orf_ids.map(_._2._2))))
      .foreach { case (node, (minimizer_size, seq_size)) => {
        pw_nii.println(node + "\t" + minimizer_size + "\t" + seq_size)
      }
      }

    pw_nii.close()
    pw_nki.close()
    pw_nki_distribution.close()
    pw_iki.close()
    pw_ni.close()
    println(timeStamp + "Successfully completed!")
  }

  def mean: Iterable[Int] => Int = values => (values.sum.toDouble / values.size).toInt
}
