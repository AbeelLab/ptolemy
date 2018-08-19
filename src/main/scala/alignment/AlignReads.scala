package alignment

/**
  * Author: Alex N. Salazar
  * Created on 15-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

import java.io.{BufferedReader, File, FileReader, PrintWriter}

import atk.Tool.progress
import utilities.FileHandling.{getFileName, getParentDirectory, openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.{NumericalUtils, SequenceUtils}
import utilities.GFAutils.ReadGFA

object AlignReads extends ReadGFA with SequenceUtils with NumericalUtils {

  case class Config(
                     canonicalQuiver: File = null,
                     reads: File = null,
                     db: File = null,
                     kmerSize: Int = 11,
                     windowSize: Int = 3,
                     proportion: Double = 0.6,
                     maxError: Double = 0.25,
                     freeMem: Int = 256,
                     minReadLength: Int = 500,
                     minHits: Int = 3,
                     verbose: Boolean = false,
                     debug: Boolean = false,
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("align-reads") {
      opt[File]('r', "reads") required() action { (x, c) =>
        c.copy(reads = x)
      } text ("Reads in FASTA-format.")
      opt[File]('c', "canonical-quiver") required() action { (x, c) =>
        c.copy(canonicalQuiver = x)
      } text ("Output directory for database to be stored.")
      opt[File]("db") required() action { (x, c) =>
        c.copy(db = x)
      } text ("Output directory for database to be stored.")
      note("\nOPTIONAL\n")
      opt[Int]("minimum-length") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Minimum read length (default is 500). Only read length of at least this size are processed.")
      opt[Int]('k', "kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Size of kmers (default is 11).")
      opt[Int]('w', "window-size") action { (x, c) =>
        c.copy(windowSize = x)
      } text ("Size of minimizer window (default is 5).")
      opt[Double]('e', "max-error") action { (x, c) =>
        c.copy(maxError = x)
      } text ("Maximum tolerable error for a read (default is 0.25). This influences the expected jaccard index " +
        "between a read and node.")
      opt[Double]('d', "distance-proportion") action { (x, c) =>
        c.copy(proportion = x)
      } text ("ORF proportion for maximum distance allowed between two minimizers during clustering. This is " +
        "based on the the size of an ORF. For example, given distance proportion, d, and size of an ORF, l, then the " +
        "max distance is: d*l.")
      opt[Int]("min-hits") action { (x, c) =>
        c.copy(minHits = x)
      } text ("Minimum number of minimizer hits for candidate node alignments (default is 3).")
      opt[Int]("free-memory") action { (x, c) =>
        c.copy(freeMem = x)
      } text ("Free memory threshold for a read chunk (default is 256).")
      opt[Unit]("verbose") action { (x, c) =>
        c.copy(verbose = true)
      }
      opt[Unit]("debug") action { (x, c) =>
        c.copy(debug = true)
      }
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.db)
      verifyFile(config.canonicalQuiver)
      alignReads(config)
    }
  }

  def alignReads(config: Config): Unit = {
    val verbose = config.verbose || config.debug
    //runtime object for getting memory usage
    val runtime = Runtime.getRuntime
    //gb in mb
    val mb = 1024 * 1024

    /**
      * Function to format memory in mb given a long value
      * @return Long
      */
    def memoryOccupation: Long => Long = x => x / mb

    //construct global kmer and node index
    val (global_kmer_index, global_node_index, global_node_info_index) = {
      //get name of canonical quiver file
      val name = getFileName(config.canonicalQuiver)
      //get parent directory
      val parent_directory = getParentDirectory(config.canonicalQuiver)
      //fetch kmer index file
      val kmer_index_file = parent_directory.listFiles().find(_.getName == name + ".ki")
      assert(kmer_index_file != None, "Could not find kmer index file for canonical quiver")
      //fetch node index file
      val node_index_file = parent_directory.listFiles().find(_.getName == name + ".ni")
      assert(node_index_file != None, "Could not find node index file for canonical quiver")
      //fetch node info index file
      val node_info_index_file = parent_directory.listFiles().find(_.getName == name + ".nii")
      assert(node_index_file != None, "Could not find node info index file for canonical quiver")
      println(timeStamp + "Loading global kmer index")
      //construct global kmer index in the format of Map[kmer, Set[Node id]]:
      // iterate through each line in the file
      val gki = openFileWithIterator(kmer_index_file.get).foldLeft(Map[Int, Set[Int]]())((kmap, line) => {
        val line_split = line.split("\t")
        //first element is assumed to be the kmer hashvalue>
        val kmer = line_split.head.toInt
        //sanity check
        assert(kmap.get(kmer) == None, "Multiple occurances of kmer hash value: " + kmer + " line: " + line)
        //construct kmer's node set, every entry is a node id
        kmap + (kmer -> line_split(1).split(",").foldLeft(Set[Int]())((node_set, node_id) => node_set + (node_id.toInt)))
      })
      println(timeStamp + "Loading global node index")
      //construct global node map in the format of Map[node id, MinKmer]:
      // iterate through each line in the file
      val gni = openFileWithIterator(node_index_file.get).foldLeft(Map[Int, Seq[MinKmer]]())((nmap, line) => {
        //split line
        val tmp = line.split(":")
        //first element assumed to be node id
        val node_id = tmp.head.toInt
        //starting positions
        val minimizer = tmp(1).split(",").map(x => {
          val tmp = x.split("_")
          new MinKmer(tmp.head.toInt, tmp(1).toInt, tmp(2).toInt)
        }).toSeq
        assert(nmap.get(node_id) == None, "Multiple occurrances of node " + node_id + ": " + line)
        nmap + (node_id -> minimizer)
      })
      //construct global node info index in the format of Map[node id, (avg minimizer set size, avg seq size)]
      val gnii = openFileWithIterator(node_info_index_file.get).foldLeft(Map[Int, (Int, Int)]())((nimap, line) => {
        //split line
        val tmp = line.split("\t")
        //get node id, average size of minimizer set, and average sequence size
        val (node_id, minimizer_set_size, seq_size) = (tmp(0).toInt, tmp(1).toInt, tmp(2).toInt)
        assert(nimap.get(node_id) == None, "Multiple occurrances of node " + node_id + ": " + line)
        nimap + (node_id -> (minimizer_set_size, seq_size))
      })
      (gki, gni, gnii)
    }

    /**
      * The expected jaccard index based on the maximum tolerable error of a read and kmer size
      */
    val expected_ji = 1 / (2 * Math.pow(Math.E, (config.maxError * config.kmerSize)) - 1)
    println(timeStamp + "Expected Jaccard Index using maximum tolerable error rate of " + (config.maxError * 100) +
      "% " + "is " + (expected_ji * 100) + "%")


    /**
      * Method to perform single-linkage clustering of minimizer hits. Returns clusters of minimizers hits as a
      * 2-tuple of (minimizer, forward position).
      *
      * @param minimizers
      * @param seq_size Size of read
      * @param max_dist
      * @return Seq[Set[(MinKmer,Int)]]
      **/
    def clusterMinimizers(minimizers: List[MinKmer], seq_size: Int, max_dist: Int): Seq[Set[(MinKmer, Int)]] = {
      /**
        * Tail-recursive implementation of clustering procedure
        *
        * @param remaining
        * @param last Last position observed in forward strand
        * @param acc
        * @param clusters
        * @return
        */
      def _clusterMinimizers(remaining: List[(MinKmer, Int)], last: Int, acc: Set[(MinKmer, Int)],
                             clusters: Seq[Set[(MinKmer, Int)]]): Seq[Set[(MinKmer, Int)]] = {
        remaining match {
          //clustered all minimizers, add last cluster if it exists
          case Nil => if (acc.isEmpty) clusters else clusters.:+(acc)
          //add to current cluster or start a new one
          case head :: tail => {
            //distance between current and last observed minimizer
            val distance = head._2 - last
            //check if distance is small enough
            if (distance <= max_dist) _clusterMinimizers(tail, head._2, acc + head, clusters)
            //start a new cluster
            else _clusterMinimizers(tail, head._2, Set(head), clusters.:+(acc))
          }
        }
      }

      //create 2-tuple of (minimizer, position in forward strand) and sort in ascending order
      val forward_map_minimizers = minimizers.map(x => {
        //compute forward position of current minimizer
        val forward_position = {
          if (x.orientation == 0) x.position else reverseComplementCoord(x.position, seq_size, config.kmerSize)
        }
        (x, forward_position)
      }).sortBy(_._2)
      if (config.debug) println(timeStamp + "----Clustering the following minimizers: " + forward_map_minimizers)
      _clusterMinimizers(forward_map_minimizers, minimizers.head.position, Set(), Seq())
    }

    //baseline memory
    val baseline_memory = memoryOccupation(runtime.totalMemory()) - memoryOccupation(runtime.freeMemory())
    println(timeStamp + "Baseline memory: " + baseline_memory + " mb")
    /**
      * Case class for FASTA-sequence entry
      * @param name
      * @param sequence
      */
    case class FastaEntry(name: String, sequence: Array[Byte])

    /**
      * Method to load set a chunk of reads up to some memory usage.
      * @param iterator
      * @param acc_iterator
      * @param acc_name
      * @param acc_sequence
      * @return List[FastaEntry]
      */
    def nextEntries(iterator: Iterator[String],
                    acc_iterator: List[FastaEntry],
                    acc_name: Option[String],
                    acc_sequence: StringBuilder): List[FastaEntry] = {
      //empty iterator or filled up available memory
      if (!iterator.hasNext || memoryOccupation(runtime.freeMemory()) < config.freeMem) {
        //no additional reads or read is below minimum read length requirement
        if (acc_name == None || acc_sequence.length < config.minReadLength) acc_iterator
          //add to chunk
        else acc_iterator.:+(new FastaEntry(acc_name.get, acc_sequence.toArray.map(encode(_))))
      }
      else {
        //get the next line
        val line = iterator.next()
        //for when the FASTA entry is empty, line is expected to be a new FASTA entry
        if (acc_name == None) {
          //sanity checks
          assert(line.startsWith(">"), "Expected start of a FASTA sequence: " + line)
          assert(acc_sequence.isEmpty, "Expected empty sequence accumulator for " + line)
          //add as start of new FASTA entry
          nextEntries(iterator, acc_iterator, Option(line.substring(1)), acc_sequence)
        }
        //for when there is an existing FASTA-entry and the next line is a new FASTA-entry
        else if (line.startsWith(">")) {
          //sanity check
          assert(acc_sequence != None, "Expected FASTA sequence for entry " + acc_name.get + ", instead" +
            " found start of new entry: " + line)
          //only add if fasta entry meets minimum read length
          if (acc_sequence.length < config.minReadLength)
            nextEntries(iterator, acc_iterator, Option(line.substring(1)), new StringBuilder)
          else {
            //create a new FastaEntry object
            val fasta_entry = new FastaEntry(acc_name.get, acc_sequence.toArray.map(encode(_)))
            //add accordingly
            nextEntries(iterator, acc_iterator.:+(fasta_entry), Option(line.substring(1)), new StringBuilder)
          }
        }
        //for when there is an existing FASTA entry and the next line is a sequence line
        else nextEntries(iterator, acc_iterator, acc_name, acc_sequence.append(line))
      }
    }
    //load reads as iterator
    val read_iterator = openFileWithIterator(config.reads)
    println(timeStamp + "Aligning reads:")
    //align reads:
    // iterate through each read
    while (read_iterator.hasNext) {
      //load next chunk of reads
      val current = nextEntries(read_iterator, List(), None, new StringBuilder)
      println(timeStamp + "--Loaded " + current.size + " reads occupying " +
        ((memoryOccupation(runtime.totalMemory()) - memoryOccupation(runtime.freeMemory())) - baseline_memory) + " mb")
      //iterate through each read in parellel and align to graph
      current.par.foreach { current_read => {
        progress(1000)
        //curent read
        //read name, seq size, and minimizers
        val (read_name, read_size, all_read_minimizers) = {
          val tmp = (current_read.name.split("\\s+").head, current_read.sequence.length,
            getMinimizers(config.kmerSize, config.windowSize, current_read.sequence))
          //the same hashvalue can have multiple minimizers
          (tmp._1, tmp._2, tmp._3.foldLeft(Map[Int, Seq[MinKmer]]())((map, minimizer) => {
            //get current minimizers for current hashvalue
            val current = map.getOrElse(minimizer.hashvalue, Seq[MinKmer]())
            map + (minimizer.hashvalue -> (current.:+(minimizer)))
          }))
        }

        if (verbose) println(timeStamp + "Processing read " + read_name + " of " + read_size + " nt")
        //get overlapping nodes and minimizers hits
        val node2kmers = {
          all_read_minimizers.foldLeft(Map[Int, Set[Int]]())((_node2kmers, minimizer) => {
            //fetch nodes that contain the current kmer
            val fetch = global_kmer_index.get(minimizer._1)
            //no nodes contain current kmer
            if (fetch == None) _node2kmers
            //add nodes to sequence of nodes
            else {
              fetch.get.foldLeft(_node2kmers)((local_node2kmers, node_id) => {
                val current = local_node2kmers.getOrElse(node_id, Set[Int]())
                local_node2kmers + (node_id -> (current + (minimizer._1)))
              })
            }
          })
            //retain nodes with enough minimizer hits
            .filter(_._2.size >= config.minHits)
        }
        //get aligned nodes in order in which it appears in the read
        val aligned_nodes = node2kmers.foldLeft(Seq[(Int, Int)]()) { case (found_nodes, (node, kmers)) => {
          if (verbose) println(timeStamp + "--Processing node " + node + " with " + kmers.size + " kmers")
          /**
            * //get all minimizers for current node
            * val all_node_minimizers = {
            * //iterate through each minimizer in the node and make a map: Node id -> Seq[Minimizers]
            * global_node_index(node).foldLeft(Map[Int, Seq[MinKmer]]())((nm_map, minimizer) => {
            * //there can be multiple minimizers with the same hashvalue
            * val current = nm_map.getOrElse(minimizer.hashvalue, Seq[MinKmer]())
            * nm_map + (minimizer.hashvalue -> (current.:+(minimizer)))
            * })
            * }
            */
          //get average minimizer set size and seq size for current node
          val (avg_ms, avg_seq_size) = global_node_info_index(node)
          //compute maximum distance threshold
          val max_dist = (avg_seq_size * config.proportion).toInt
          //get the most dense cluster that are epsilon-based chained
          val (max_dense_cluster, total_clusters) = {
            //cluster minimizers on the read
            val tmp = clusterMinimizers(kmers.toSeq.map(all_read_minimizers(_)).flatten.toList, read_size, max_dist)
            if (verbose) println(timeStamp + "----Found " + tmp.size + " clusters")
            if (config.debug) println(timeStamp + "----Formed the following clusters: " + tmp)
            (if (tmp.isEmpty) Seq() else tmp.maxBy(_.size), tmp.size)
          }
          if (verbose)
            println(timeStamp + "----Using max dense cluster of " + max_dense_cluster.size + " " + "minimizers")
          //compute jaccard index
          val jaccard_index = max_dense_cluster.size.toDouble / avg_ms
          if (verbose)
            println(timeStamp + "----Jaccard index of " + jaccard_index + " based on mean minimizer set size of " + avg_ms)
          //jaccard index of node is below expected, no alignment
          if (jaccard_index < expected_ji) found_nodes
          //jaccard index is at least of the expected score, consider alignment
          else {
            //compute median position
            val median_position = computeMedian(max_dense_cluster.map(_._2))
            if (verbose)
              println(timeStamp + "----Node is alignable. Using median coordinate position of: " + median_position)
            found_nodes.:+((node, median_position))
          }
        }
        }.sortBy(_._2)

        if (verbose) {
          val message = if (aligned_nodes.isEmpty) "--No alignments found" else "--Found " + aligned_nodes.size + " node" +
            " alignments in the following order: " + aligned_nodes
          println(timeStamp + message)
        }
      }
      }
    }
  }
}
