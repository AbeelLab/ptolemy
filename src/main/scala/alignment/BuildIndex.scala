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
import utilities.SequenceUtils.encode

object BuildIndex extends ReadGFA with PtolemyDB with MMethods{

  case class Config(
                     canonicalQuiver: File = null,
                     db: File = null,
                     kmerSize: Int = 11,
                     windowSize: Int = 3,
                     maxNodeFreq: Double = 0.01,
                     verbose: Boolean = false
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("build-index") {
      opt[File]('c', "canonical-quiver") required() action { (x, c) =>
        c.copy(canonicalQuiver = x)
      } text ("Path to canonical quiver in GFA-format.")
      opt[File]("db") required() action { (x, c) =>
        c.copy(db = x)
      } text ("Directory path of database (e.g. output directory of 'extract' module).")
      note("\nOPTIONAL\n")
      opt[Int]('k', "kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Size of kmers (default is 11).")
      opt[Int]('w', "window-size") action { (x, c) =>
        c.copy(windowSize = x)
      } text ("Size of minimizer window (default is 5).")
      opt[Double]("max-frequency") action { (x, c) =>
        c.copy(maxNodeFreq = x)
      } text ("Ignore kmers that appear in this .")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      }
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.db)
      verifyFile(config.canonicalQuiver)
      buildIndex(config)
    }
  }

  def buildIndex(config: Config): Unit = {
    println(timeStamp + "Loading canonical quiver")
    //load all nodes in the canonical quiver
    val gfa_nodes = loadNodesGFA(config.canonicalQuiver)
    println(timeStamp + "Found " + gfa_nodes.size + " nodes in the canonical quiver")
    println(timeStamp + "Loading ORF ID to node ID schema")
    //create map from orf ID to node IDs
    val orfid2nodeid = loadORFid2Nodeid(config.db)
    //sanity check
    gfa_nodes.foreach(x => assert(!orfid2nodeid.get(x).isEmpty, "Could not find mapping for " + x))

    /**
      * Method create a minimizer hashtable and node info hashtable given fasta file and empty/existing hashtables.
      * The minimizer hashtable schema: kmer -> Map(node ID -> Set(position, orientation))
      * The node info hashtable schema: node id -> Map(orf id -> (minimizer set size, seq size))
      *
      *
      * @param fasta_file File object of fasta file
      * @param kmer_index Empty or existing hashtable
      * @param nodeinfo_map Empty or existing hashtable
      * @return (Map[Int, Map[Int, Set[(Int,Int)], Map[Int, Map[Int, (Int,Int)])
      **/
    def fetchMinimizers(fasta_file: File,
                        kmer_index: Map[Int, Map[Int, Set[(Int,Int)]]],
                        nodeinfo_map: Map[Int, Map[Int, (Int,Int)]]
                       ): (Map[Int, Map[Int, Set[(Int,Int)]]], Map[Int, Map[Int, (Int,Int)]]) = {
      //create fasta iterator
      val fasta_iterator = new FastaIterator(fasta_file)
      /**
        * Tail-recusive method to iterate through fasta entries and update hashtable
        *
        * @param _kmer_index Hashtable to update
        * @return Map[Int, Map[Int,Int]
        **/
      def _fetchMinimizers(_kmer_index: Map[Int, Map[Int, Set[(Int,Int)]]],
                          _nodeinfo_map: Map[Int, Map[Int, (Int,Int)]]
                          ): (Map[Int, Map[Int, Set[(Int,Int)]]], Map[Int, Map[Int, (Int,Int)]]) = {
        if (!fasta_iterator.hasNext) (_kmer_index, _nodeinfo_map)
        else {
          //get current orf
          val current_orf = fasta_iterator.next()
          //get orf id and sequence
          val (orf_id, orf_seq) = (current_orf.getDescription.substring(1).toInt, current_orf.getSequence)
          //get corresponding node id
          val node_id = orfid2nodeid(orf_id)
          //get minimizers
          val minimizers = getMinimizers(config.kmerSize, config.windowSize, orf_seq.toCharArray.map(encode(_)))
          //update node info map
          val updated_nodeinfo = {
            val current = _nodeinfo_map.getOrElse(node_id, Map[Int, (Int,Int)]())
            assert(current.get(orf_id) == None, "Multiple occurrences of ORF sequence: " + orf_id)
            _nodeinfo_map + (node_id -> (current + (orf_id -> (minimizers.size, orf_seq.size))))
          }
          //update kmer index
          _fetchMinimizers(minimizers.foldLeft(_kmer_index)((local_kmer_index, minimizer) => {
            //get kmer's map of node ids to starting positions and update accordingly
            val current = {
              //get all nodes where this kmer appears in
              val node_map = local_kmer_index.getOrElse(minimizer.hashvalue, Map[Int, Set[(Int,Int)]]())
              //get all starting positions of current node
              val node_seq = node_map.getOrElse(node_id, Set[(Int,Int)]())
              //update starting positions of current node
              node_map + (node_id -> (node_seq + ((minimizer.position, minimizer.orientation))))
            }
            //update kmer index
            local_kmer_index + (minimizer.hashvalue -> current)
          }),updated_nodeinfo)
        }
      }
      _fetchMinimizers(kmer_index, nodeinfo_map)
    }

    println(timeStamp + "Creating global kmer index")
    //create a global kmer index: hashtable with kmer as keys and a map of node id -> start position as values
    // as well as global node info index: hashtable with node id as keys and a map of
    // orf id -> (avg minimizer setsize, avg sequence size)
    val (global_kmer_index, global_nodeinfo) = {
      //fetch all ORF sequences files
      config.db.listFiles().filter(_.getName.endsWith("fasta"))
        //iterate through each one and update kmer index (Map[kmer -> Map[Node, starting positions]])
        .foldLeft((Map[Int, Map[Int, Set[(Int,Int)]]](), Map[Int, Map[Int, (Int,Int)]]())) {
        case ((kmer_index, nodeinfo), _fasta) => fetchMinimizers(_fasta, kmer_index, nodeinfo) }
    }

    println(timeStamp + "Creating global node index")
    //using the global kmer index, create hashtable with node ids as keys and ordered sequence of starting positions
    val global_node_index = {
      //iterate through each kmer and it's map of node id's, ultimately updating the node index
      global_kmer_index.foldLeft(Map[Int, Set[(Int,Int,Int)]]()) { case (node_index, (kmer, node_map)) => {
        //iterate through each node in the node map and update the global starting positions
        node_map.foldLeft(node_index) { case (local_node_index, (node_id, starting_positions)) => {
          //get current starting positions
          val current_starting_positions = local_node_index.getOrElse(node_id, Set[(Int,Int,Int)]())
          //update starting positions
          local_node_index + (node_id -> (current_starting_positions ++ starting_positions.map(x => (kmer,x._1,x._2))))
        }}
      }}
      //sort starting positions in ascending order
        .mapValues(_.toSeq.sortBy(x => (x._2, x._1)))
    }
    println(timeStamp + "Writing to disk")
    //fetch name of canonical quiver file
    val name = getFileName(config.canonicalQuiver)
    //fetch parent directory
    val parent_directory = getParentDirectory(config.canonicalQuiver)
    //create output files for db indeces
    val pw_ki = new PrintWriter(parent_directory + "/" + name + ".ki")
    val pw_ni = new PrintWriter(parent_directory + "/" + name + ".ni")
    val pw_nii = new PrintWriter(parent_directory + "/" + name + ".nii")
    //output kmer index
    // this is different the data structure constructed above, it's structured: <kmer hash>\t<node id>,<node id>...
    global_kmer_index.foreach { case (kmer, node_map) =>
      pw_ki.println(kmer + "\t" + node_map.foldLeft(Set[Int]())((b,a) => b + (a._1)).mkString(","))
    }
    //output node index
    global_node_index.foreach { case (node, starting_positions) =>
      pw_ni.println(node + ":" + starting_positions.map(x => (x._1 + "_" + x._2 + "_" + x._3)).mkString(","))
  }
    //summarize node info as: node id -> (average minimizer set size, average sequence size)
    global_nodeinfo.mapValues(orf_ids => (mean(orf_ids.map(_._2._1)), mean(orf_ids.map(_._2._2))))
      .foreach{case (node, (minimizer_size, seq_size)) => {
        pw_nii.println(node + "\t" + minimizer_size + "\t" + seq_size)
      }}

    pw_nii.close()
    pw_ki.close()
    pw_ni.close()
    println(timeStamp + "Successfully completed!")
  }

  def mean: Iterable[Int] => Int = values => (values.sum.toDouble / values.size).toInt
}
