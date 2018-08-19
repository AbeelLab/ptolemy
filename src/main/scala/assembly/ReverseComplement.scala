package assembly

/**
  * Author: Alex N. Salazar
  * Created on 14-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

import java.io.{File, PrintWriter}
import utilities.FileHandling._
import utilities.SequenceUtils


object ReverseComplement {

  case class Config(
                     fasta: File = null,
                     gff: File = null,
                     outputDir: File = null,
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("syntenic-anchors") {
      opt[File]('a', "assembly") required() action { (x, c) =>
        c.copy(fasta = x)
      } text ("Path to assembly is FASTA format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory.")
      note("\nOPTIONAL\n")
      opt[File]("gff") action { (x, c) =>
        c.copy(gff = x)
      } text ("GFF-formatted file. Note: only changes coordinates.")

    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyFile(config.fasta)
      verifyDirectory(config.outputDir)
      makeReverseComplement(config)
    }
  }

  def makeReverseComplement(config: Config): Unit = {
    println(timeStamp + "Creating reverse complement of assembly")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" +
      config.fasta.getName.substring(0, config.fasta.getName.lastIndexOf(".")) + ".reverse_complement.fasta")
    val empty: Option[(String, String)] = None
    //iterate through fasta file construct map of assembly size and output reverse complements
    val assembly_size = {
      val tmp = openFileWithIterator(config.fasta).foldLeft((Map[String, Int](), empty)) {
        case ((assembly_size, acc), line) => {
          //very first line of FASTA file
          if (acc == None) {
            assert(line.startsWith(">"), "Expected start of FASTA entry: " + line)
            //output line
            pw.println(line)
            //return map, add fasta entry and add assembly size
            (assembly_size, Some(line.substring(1).split("\\s+").head, ""))
          }
          //start of new fasta entry
          else if (line.startsWith(">")) {
            //output reverse entry of fasta entry
            acc.get._2.reverseIterator.foreach(nt => pw.print(reverseComplement(nt)))
            pw.println
            pw.println(line)
            //add assembly size of previous entry, start new acc
            (assembly_size + (acc.get._1 -> acc.get._2.size), Some(line.substring(1).split("\\s+").head, ""))
          }
          //sequence line from continuing fasta entry, append line to sequence
          else (assembly_size, Some(acc.get._1, acc.get._2 + line))
        }
      }
      //add remaining fasta entry if existing
      if (tmp._2 == None) tmp._1
      else {
        //output reverse entry of fasta entry
        tmp._2.get._2.reverseIterator.foreach(nt => pw.print(reverseComplement(nt)))
        pw.println
        tmp._1 + (tmp._2.get._1 -> tmp._2.get._2.size)
      }
    }
    pw.close()
    //update GFF file if provided
    if (config.gff != null) {
      println(timeStamp + "Creating reverse complement of GFF file")
      //output file
      val pw2 = new PrintWriter(config.outputDir + "/" +
        config.gff.getName.substring(0, config.gff.getName.lastIndexOf(".")) + ".reverse_complement.gff")
      //iterate through GFF file and return reverse complemented coordinates
      val (header, updated_gff) =
      openFileWithIterator(config.gff).foldLeft(Seq[String](), Seq[Array[String]]())((_updated_gff, line) => {
        if (line.startsWith("#")) (_updated_gff._1.:+(line), _updated_gff._2)
        else {
          val tmp = line.split("\t").zipWithIndex
          //line with reverse complemented coordinates
          val array = tmp.map(x => {
            //if not the start or end columns, move on
            if (x._2 != 3 && x._2 != 4) x._1
            //get reverse complement coordinate
            else reverseComplementCoord(x._1.toInt, assembly_size(tmp.head._1)).toString
          })
          (_updated_gff._1, _updated_gff._2.:+(array.updated(3, array(4)).updated(4, array(3))))
        }
      })
      header.foreach(pw2.println)
      updated_gff.sortBy(x => (x(3).toInt, x(4).toInt)).foreach(x => pw2.println(x.mkString("\t")))
      pw2.close
    }
    println(timeStamp + "Successfully completed!")
  }

  /**
    * Function to compute the corresponding reverse complement coordinate for a given coordinate
    *
    * @return Int
    */
  def reverseComplementCoord: (Int, Int) => Int = (position, seq_size) => seq_size - (position - 1)

  /**
    * IUPAC-supported reverse complement function
    * @return Char
    */
  def reverseComplement: Char => Char = nt => {
    nt match {
      case 'A' => 'T'
      case 'a' => 't'
      case 'T' => 'A'
      case 't' => 'a'
      case 'G' => 'C'
      case 'g' => 'c'
      case 'C' => 'G'
      case 'c' => 'g'
      case 'Y' => 'R'
      case 'y' => 'r'
      case 'R' => 'Y'
      case 'r' => 'y'
      case 'S' => 'S'
      case 's' => 's'
      case 'W' => 'W'
      case 'w' => 'w'
      case 'K' => 'M'
      case 'k' => 'm'
      case 'M' => 'K'
      case 'm' => 'k'
      case 'B' => 'V'
      case 'b' => 'v'
      case 'D' => 'H'
      case 'd' => 'h'
      case 'H' => 'D'
      case 'h' => 'd'
      case 'V' => 'B'
      case 'v' => 'b'
      case 'N' => 'N'
      case 'n' => 'n'
      case _ => {
        assert(false, "Could not recognize the following nucleotide: " + nt);
        'Z'
      }
    }
  }

}
