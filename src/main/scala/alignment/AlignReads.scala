package alignment

/**
  * Author: Alex N. Salazar
  * Created on 15-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

import java.io.{File, PrintWriter}

import atk.Tool.progress
import utilities.FileHandling._
import utilities.NumericalUtils
import utilities.SequenceUtils.{FastaUtils}
import utilities.DatabaseUtils.GraphIndex
import utilities.MinimizerUtils.{Minimizer, MMethods}
import utilities.GFAutils.ReadGFA

object AlignReads extends ReadGFA with NumericalUtils with GraphIndex with FastaUtils with MMethods {

  case class Config(
                     canonicalQuiver: File = null,
                     reads: File = null,
                     db: File = null,
                     kmerSize: Int = 11,
                     windowSize: Int = 3,
                     proportion: Double = 0.6,
                     maxError: Double = 0.25,
                     chunkSize: Int = 20000,
                     minReadLength: Int = 500,
                     minHits: Int = 3,
                     verbose: Boolean = false,
                     prefix: String = null,
                     outputDir: File = null,
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
      } text ("Directory path of database (e.g. output directory of 'extract' module).")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory for database to be stored.")
      opt[String]('p', "prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
      note("\nOPTIONAL\n")
      opt[Int]("minimum-length") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Minimum read length (default is 500). Only read length of at least this size are processed.")
      opt[Int]('k', "kmer-size") action { (x, c) =>
        c.copy(kmerSize = x)
      } text ("Size of kmers (default is 11).")
      opt[Int]('w', "window-size") action { (x, c) =>
        c.copy(windowSize = x)
      } text ("Size of minimizer window (default is 3).")
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
      opt[Int]("chunk-size") action { (x, c) =>
        c.copy(chunkSize = x)
      } text ("Number of reads to load at a time(default is 30000).")
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

  /**
    * The expected jaccard index based on the maximum tolerable error of a read and kmer size
    *
    * @return Double
    */
  def expectedJI: (Double, Int) => Double = (error, kmer_size) => 1 / (2 * Math.pow(Math.E, (error * kmer_size)) - 1)

  def updateCanonicalQuiver(output_directory: File, output_name: String, cq: File): Unit = {
    val tmp_output = new File(output_directory + "/." + output_name + ".tmp")
    println(timeStamp + "Computing edge coverage")
    //iterate through alignments and keep track of edges and their coverage as well as total reads and unmapped
    val (edge_coverage, node_coverage, total_reads, mapped) =
    //itereate each line with edge coverage, node coverage, total reads ,total mapped
      openFileWithIterator(tmp_output).foldLeft((Map[(Int, Int), Int](), Map[Int, Int](), 0, 0)) {
        case ((ecoverage, ncoverage, n, m), line) => {
          //parse alignment line
          val (read_id, alignment) = parseAlignment(line)
          //if unmapped or maps to only one node, do not update edge coverage
          ((if (alignment.size < 2) ecoverage
          //else updated coverage of observed edges as well as total number of reads
          else alignment.sliding(2).foldLeft(ecoverage)((local_ecoverage, _edge) => {
            val edge = (_edge(0), _edge(1))
            val current_coverage = local_ecoverage.getOrElse(edge, 0)
            local_ecoverage + (edge -> (current_coverage + 1))
          }))
            //incremenet node coverages
            , (alignment.foldLeft(ncoverage)((local_ncoverage, node) => {
            val current_coverage = local_ncoverage.getOrElse(node, 0)
            local_ncoverage + (node -> (current_coverage + 1))
          }))
            //increment total reads and total mapped
            , n + 1, m + 1)
        }
      }
    println(timeStamp + "Found " + total_reads + " total reads, of which " + (mapped.toDouble / total_reads) * 100 + "% are mapped")
    println(timeStamp + "Writing coverage to disk")
    val pw = new PrintWriter(output_directory + "/" + output_name)
    //load canonical quiver and add coverage information while retaining new edges
    openFileWithIterator(cq).foldLeft(edge_coverage)((filtered_coverage, line) => {
      //get line type
      val line_type = line.split("\t").head
      line_type match {
        //edge line
        case "L" => {
          //get edge
          val edge = parseLinkLine(line)
          //get coverage
          val coverage = edge_coverage.getOrElse(edge, 0)
          pw.println(line + "\tRC:" + coverage)
          if (coverage == 0) filtered_coverage else filtered_coverage - edge
        }
        //node line
        case "S" => {
          //get node
          val node = parseSegmentLine(line)
          //get coverage
          val coverage = node_coverage.getOrElse(node, 0)
          pw.println(line + "\tRC:" + coverage)
          filtered_coverage
        }
        //Something else, move on
        case _ => {
          pw.println(line)
          filtered_coverage
        }
      }
    })
      //add new edges to output
      .foldLeft(0)((index, new_edge) => {
      if (index == 0) pw.println("#NEW EDGES FROM READ ALIGNMENTS")
      pw.println("L\t" + new_edge._1._1 + "\t+\t" + new_edge._1._2 + "\t+\t1M\tRC:" + new_edge._2)
      index + 1
    })
    //iterate through alignments once more and output to file
    pw.println("#READ ALIGNMENTS")
    openFileWithIterator(tmp_output).foreach(pw.println)
    pw.close()
    tmp_output.delete()
    println("Successfully completed!")
  }

  def alignReads(config: Config): Unit = {
    //set log level
    val verbose = config.verbose || config.debug
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
      val gki = loadGlobalKmerIndex(kmer_index_file.get)
      println(timeStamp + "Loading global node index")
      //construct global node map in the format of Map[node id, MinKmer]:
      val gni = loadGlobalNodeIndex(node_index_file.get)
      println(timeStamp + "Loading global node info index")
      //construct global node info index in the format of Map[node id, (avg minimizer set size, avg seq size)]
      val gnii = loadGlobanNodeInfoIndex(node_info_index_file.get)
      //return indeces
      (gki, gni, gnii)
    }
    //compute expected jaccard index for the given max tolerable error rate and kmer size
    val expected_ji = expectedJI(config.maxError, config.kmerSize)
    println(timeStamp + "Expected Jaccard Index using maximum tolerable error rate of " + (config.maxError * 100) +
      "% " + "is " + (expected_ji * 100) + "%")

    //create temporary file path
    val output_name = config.prefix + ".gfa"
    //create output file
    val pw = new PrintWriter(config.outputDir + "/." + output_name + ".tmp")
    //load reads as iterator
    val read_iterator = openFileWithIterator(config.reads)
    var next_line: Option[String] = None
    println(timeStamp + "Starting alignments")
    //iterate through each read and align
    while (read_iterator.hasNext) {
      //load next chunk of reads
      val (current, runoff_line) = loadFastaChunk(read_iterator, List(), config.chunkSize, config.minReadLength,
        next_line, new StringBuilder)
      //get carry over read entry
      next_line = runoff_line
      println(timeStamp + "Loaded " + current.size + " reads")
      //iterate through each read in parellel and align to graph
      current.foreach { current_read => {
        progress(1000)
        //curent read
        //read name, seq size, and minimizers
        val (read_name, read_size, all_read_minimizers) = {
          val tmp = (current_read.name.split("\\s+").head, current_read.sequence.length,
            getMinimizers(config.kmerSize, config.windowSize, current_read.sequence))
          //the same hashvalue can have multiple minimizers
          (tmp._1, tmp._2, tmp._3.foldLeft(Map[Int, Seq[Minimizer]]())((map, minimizer) => {
            //get current minimizers for current hashvalue
            val current = map.getOrElse(minimizer.hashvalue, Seq[Minimizer]())
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
        //get aligned nodes in order in which it appears in the read in 2-tuple (node, position, orientation)
        // note: orientation is simply majority count of the read-node minimizer orientation occurrance
        val aligned_nodes = node2kmers.foldLeft(Seq[(Int, (Int, Int), Int)]()) { case (found_nodes, (node, kmers)) => {
          if (verbose) println(timeStamp + "--Processing node " + node + " with " + kmers.size + " kmers")
          //get average minimizer set size and seq size for current node
          val (avg_ms, avg_seq_size) = global_node_info_index(node)
          //compute maximum distance threshold
          val max_dist = (avg_seq_size * config.proportion).toInt
          //get the most dense cluster that are epsilon-based chained
          val (max_dense_cluster, total_clusters) = {
            //cluster minimizers on the read
            val tmp = clusterMinimizers(kmers.toSeq.map(all_read_minimizers(_)).flatten.toList,
              config.kmerSize, read_size, max_dist, config.debug)
            if (verbose) println(timeStamp + "----Found " + tmp.size + " clusters")
            if (config.debug) println(timeStamp + "----Formed the following clusters: " + tmp)
            (if (tmp.isEmpty) Set[(Minimizer, Int)]() else tmp.maxBy(_.size).toSet, tmp.size)
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
            //create a map of the minimizer hashvalue -> orientation(s)
            val read_kmers2orientation = max_dense_cluster.foldLeft(Map[Int, Seq[Int]]()) {
              case (ro_map, (minimizer, position)) => {
                val current = ro_map.getOrElse(minimizer.hashvalue, Seq[Int]())
                ro_map + (minimizer.hashvalue -> (current.:+(minimizer.orientation)))
              }
            }
            //pair the orientations of a read minimizer hit to the node minimizer hit and count their occurrances
            // as well as creating a 2-tuple of (read-forward position, node-forward position)
            val (orientation, start, end) = {
              //fetch all node minimizers and iterate through them with map of orientation occurance -> count
              val tmp = global_node_index(node).foldLeft((Map[Int, Int](), Seq[(Int, Int)]())) {
                case ((zo_map, read2node_positions), minimizer) => {
                  //attempt to get the corresponding orientation of current minimizer in the read
                  val read_orientations = read_kmers2orientation.get(minimizer.hashvalue)
                  //minimizer does not exist in read, move on
                  if (read_orientations.isEmpty) (zo_map, read2node_positions)
                  else {
                    val node_forward_position = forwardPosition(minimizer, avg_seq_size, config.kmerSize)
                    //iterate through each orientation
                    (read_orientations.get.foldLeft(zo_map)((local_zo_map, read_orientation) => {
                      //determine orientation of minimizer instance
                      val orientation_instance = getOrientation(read_orientation, minimizer.orientation)
                      //get the current count of the orientation occurrance
                      val current = local_zo_map.getOrElse(orientation_instance, 0)
                      //increment accordingly
                      local_zo_map + (orientation_instance -> (current + 1))
                    }),
                      //get corresponding read minimizer(s) and add forward positions to sequence
                      all_read_minimizers(minimizer.hashvalue).foldLeft(read2node_positions)((local_read2node_pos, rm) => {
                        local_read2node_pos.:+((forwardPosition(rm, read_size, config.kmerSize), node_forward_position))
                      }))
                  }
                }
              }
              val sorted = tmp._2.sortBy(_._2)

              //get orientation of alignment, this is just to be able to output occurrance map in debug mode
              val final_orientation = tmp._1.toList.sortBy(-_._2).head._1

              if(config.debug) println(timeStamp + "----Orientation map: " + tmp._1)

              def getBest(_sorted: Seq[(Int, Int)]): (Int,Int) = {
                if (_sorted.isEmpty) {
                  assert(false, "Could not find reliable coordinates for read " + read_name + " during alignment with" +
                    " node " + node + " using max dense cluster of " + max_dense_cluster)
                  (-1,-1)
                }
                else {
                  val left_most = _sorted.head
                  val right_most = _sorted.last
                  //compute predicted star/end coordinates on the read
                  val (start, end) = {
                    //calculate according to the orientation
                    if (final_orientation == 0) (left_most._1 - left_most._2, right_most._1 + abs(avg_seq_size - right_most._2))
                    else (right_most._1 - abs(avg_seq_size - right_most._2), left_most._1 + left_most._2)
                  }
                  if(start < end) (start, end)
                  else {
                    if(config.debug)
                      println(timeStamp + "----WARNING: predicted start/end is illogical: " +(left_most, right_most))
                    val filtered_sorted = _sorted.take(sorted.size - 1).tail
                    getBest(filtered_sorted.take(filtered_sorted.size - 1).tail)
                  }
                }
              }
              val (_start, _end) = getBest(sorted)
              //sort orientation occurance by most frequent
              (final_orientation, _start, _end)
            }
            if (verbose) {
              println(timeStamp + "----Node is alignable.")
              println(timeStamp + "----Node has an average size of: " + avg_seq_size)
              println(timeStamp + "----Predicted start/end coordinates: " + (start, end))
              println(timeStamp + "----Choosing orientation: " + orientation)
            }
            found_nodes.:+((node, (start, end), orientation))
          }
        }
        }.sortBy(x => (x._2._1, x._2._2))

        //if read is unmapped
        if (aligned_nodes.isEmpty) {
          if (verbose) println(timeStamp + "--No alignments found")
          pw.println("P" + "\t" + read_name + "\t")
        }
        //read is mapped
        else {
          if (config.debug) println(timeStamp + "Found the following alignment: " + aligned_nodes)
          //attempt to determine orientation of read based on orientation of each node
          val alignment_orientation = {
            //sum total occurance of forward and reverse nodes
            val tmp = aligned_nodes.foldLeft((0, 0)) { case ((f, r), pnode) => if (pnode._3 == 0) (f + 1, r) else (f, r + 1) }
            //majority vote, unless ambiguous
            if (tmp._1 == tmp._2) '?'
            else if (tmp._1 > tmp._2) '+'
            else '-'
          }
          //get sequence of predicted node alignments
          val aligned_structure = {
            val tmp = aligned_nodes.map(_._1 + "+")
            //determine correct orientation
            if (alignment_orientation == '+') tmp else tmp.reverse
          }.mkString(",")
          if (verbose) println(timeStamp + "--Found " + aligned_nodes.size + " node alignments in the following " +
            "order: " + aligned_structure)
          pw.println("P" + "\t" + read_name + "\t" + aligned_structure)
        }
      }
      }
    }
    pw.close
    println(timeStamp + "Alignment completed")
    updateCanonicalQuiver(config.outputDir, output_name, config.canonicalQuiver)
  }
}
