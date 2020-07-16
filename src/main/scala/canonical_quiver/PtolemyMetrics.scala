package canonical_quiver

/**
  * Author: Alex N. Salazar
  * Created on 26-3-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

import java.io.{File, PrintWriter}
import java.util.Locale

import utilities.{GFFutils, MinimapUtils}
import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}

import scala.collection.immutable.HashMap

object PtolemyMetrics extends GFFutils with MinimapUtils {


  case class Config(
                     quiver: File = null,
                     outputDir: File = null,
                     bin: Double = 0.01,
                     verbose: Boolean = false,
                   )

  def main(args: Array[String]): Unit = {
    Locale.setDefault(new Locale("en", "US"))
    val parser = new scopt.OptionParser[Config]("ptolemy-metrics") {
      opt[File]('c', "canonical-quiver") required() action {
        (x, c) => c.copy(quiver = x)
      } text ("Canonical quiver representation.")
      opt[File]('o', "output-directory") required() action {
        (x, c) => c.copy(outputDir = x)
      } text ("Output directory.")
      note("\nOPTIONAL\n")
      opt[Unit]("verbose") action {
        (x, c) => c.copy(verbose = true)
      }

    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.quiver)
      ptolemyMetrics(config)
    }
  }

  def ptolemyMetrics(config: Config): Unit = {
    println(timeStamp + "Processing GFA file")
    //open canonical quiver in GFA format
    val (sequence2nodes, genome2sequence, node2labels) = openFileWithIterator(config.quiver)
      .foldLeft((HashMap[String, Seq[Int]](), HashMap[String, Set[String]](), HashMap[Int, Set[String]]())) {
        case ((_sequence2nodes, _genome2sequence, _node2labels), line) => {
          //move on if it's not a path or genome line
          if (!line.startsWith("P") && !line.startsWith("G")) (_sequence2nodes, _genome2sequence, _node2labels)
          else {
            //process path line
            if (line.startsWith("P")) {
              //get sequence and path
              val (sequence, path) = getPathLine(line)
              //assert sequence has not been seen before
              assert(_sequence2nodes.get(sequence) == None)
              (_sequence2nodes + (sequence -> path), _genome2sequence,
                //iterate through each path and update nodes and their corresponding labels
                path.foldLeft(_node2labels)((local_node2labels, node) => {
                  val get = local_node2labels.getOrElse(node, Set[String]())
                  local_node2labels + (node -> (get + sequence))
                }))
            }
            //porcess genome line
            else {
              val (genome, sequences) = getGenomeLine(line)
              assert(_genome2sequence.get(genome) == None)
              (_sequence2nodes, _genome2sequence + (genome -> sequences), _node2labels)
            }
          }
        }
      }
    //reverse the hashmap
    val sequence2genome = genome2sequence.toList.foldLeft(HashMap[String, String]()){case (map, (genome,sequences)) => {
      sequences.toList.foldLeft(map)((local_map, sequence) => {
        local_map + (sequence -> genome)
      })
    }}
    //create output file for total number of native genes
    val pw_average_genes = new PrintWriter(config.outputDir + "/native_genes_distribution.txt")
    //output header
    pw_average_genes.println("Genome\tGenes")
    pw_average_genes.println("Quiver\t" + node2labels.size)
    //create output file for relative location of structural variants
    val pw_location_sv = new PrintWriter(config.outputDir + "/relative_location_of_svs.txt")
    pw_location_sv.println("Relative_location\tLabels")
    //group by genome and calculate total number of genes per genome
    sequence2nodes.toList.groupBy(x =>sequence2genome(x._1)).foreach(genome => {
      //output total number of genes in a genome
      pw_average_genes.println(Seq(genome._1, genome._2.map(_._2.size).sum).mkString("\t"))
      //iterate through each sequence, output
      genome._2.foreach{ case (sequence_id, path) => {
        //get total nodes in path
        val total_nodes = path.size
        //iterate through each node, get relative location
        path.zipWithIndex.foldLeft(config.bin){case(current_bin,(node, index)) =>{
          val relative_location = normalizeIndex(index, total_nodes)
          if(relative_location <= current_bin) {
            pw_location_sv.println(current_bin + "\t" + node2labels(node).size)
            current_bin
          } else {
            val new_bin = "%.2f".format(current_bin + config.bin).toDouble
            pw_location_sv.println(new_bin + "\t" + node2labels(node).size)
            new_bin
          }
        }}}
      }
    })
    pw_location_sv.close
    val pw_location_sv_summarized = new PrintWriter(config.outputDir + "/relative_location_of_svs.summarized.txt")
    pw_location_sv_summarized.println("Bin\tMean\tSd")
    openFileWithIterator(new File(config.outputDir + "/relative_location_of_svs.txt")).toList.tail.map(_.split("\t"))
        .groupBy(_.head).mapValues(x => {
            val values = x.map(_(1).toInt)
            val mean = values.sum / values.size.toDouble
            val sd = Math.sqrt(values.map(x => Math.pow(x - mean, 2)).sum / (values.size - 1))
            (mean,sd)
    }).toList.sortBy(_._1).foreach(x => pw_location_sv_summarized.println(x._1 + "\t" + x._2._1 + "\t" + x._2._2))
    pw_location_sv_summarized.close
    pw_average_genes.close
    //create output file for syntenic anchor distribution
    val pw = new PrintWriter(config.outputDir + "/syntenic_anchor_size_distribution.txt")
    //output header
    pw.println("Node\tLabels")
    node2labels.foreach(x => pw.println(Seq(x._1, x._2.size).mkString("\t")))
    pw.close
    println(timeStamp + "Successfully completed!")
  }

  /**
    * Function to parse a path line from GFA file
    *
    * @return (String, Seq[Int])
    */
  def getPathLine: String => (String, Seq[Int]) = line => {
    val split = line.split("\t")
    (split(1), split(2).split("""\W+""").map(_.toInt).toSeq)
  }

  /**
    * Function to parse a genome line from GFA file
    *
    * @return (String, Seq[Int])
    */
  def getGenomeLine: String => (String, Set[String]) = line => {
    val split = line.split("\t")
    (split(1), split(2).split(",").toSet)
  }

  def normalizeIndex: (Int,Int) => Double = (index, size) => index.toDouble / size

}
