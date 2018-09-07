package utilities

import utilities.FileHandling.timeStamp
import utilities.NumericalUtils.{max, min}
import scala.annotation.tailrec
import scala.collection.immutable.HashMap

/**
  * Author: Alex N. Salazar
  * Created on 5-12-2017
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

trait GFFutils {

  /**
    * Case class for storing a GFF gene annotation as an object
    *
    * @param chrm
    * @param sample
    * @param feature
    * @param start
    * @param end
    * @param orientation
    * @param name
    * @param id
    */
  case class GFFLine(chrm: String, sample: String, feature: String,
                     start: Int, end: Int, orientation: String, name: Option[String], id: String) {
    //get parent attribute
    lazy val parent = if (feature == "gene") None else Option(id.replace("cds", "gene"))
    //get original attribute
    lazy val attribute =
      "ID=" + id +
        (if (parent != None) (";Parent=" + parent.get) else "") +
        (if (parent == None) ";Name=" + name.get else "")

    //override toString to return original GFFLine
    override def toString = {
      Seq(chrm, sample, feature, start, end, ".", orientation, ".", attribute).mkString("\t")
    }
  }

  /**
    * Function to parse gff line and return GFFLine object
    *
    * @return GFFLine object
    */
  def toGFFLine: (String, Boolean) => Option[GFFLine] = (line, showWarnings) => {
    val tmp = line.split("\t")
    if(tmp.size < 9) {
      //log message
      if (showWarnings) println(timeStamp + "----WARNING: found a non-GFF-formatted line. Skipping: " + line)
      None
    }
    else Option(new GFFLine(tmp(0), tmp(1), tmp(2), tmp(3).toInt, tmp(4).toInt,
      tmp(6), getGeneName(tmp(8)), getGeneID(tmp(8))))
  }

  /** Harcoded string in the GFF file to identify unique gene ID. Used in the method, 'getGeneID' */
  private val gene_id_string = "ID="

  /**
    * Method to parse out the gene name in the attributes column of a gff file
    */
  private def getGeneName(attribute: String): Option[String] = {
    val annotation = attribute.split(";")
    val locus_tag = annotation.find(_.contains("locus_tag="))
    val gene_name = annotation.find(_.contains("Name="))
    if(locus_tag.isEmpty){
      if(gene_name.isEmpty) None
      else Option(gene_name.get.substring(5))
    }
    else {
      val lt = locus_tag.get.substring(10)
      if(gene_name.isEmpty) Option(lt)
      else Option(lt + "_" + gene_name.get.substring(5))
    }
  }

  /**
    * Method to extract the local gene id
    */
  private def getGeneID(x: String): String = {
    val tmp = x.split(";").find(_.contains(gene_id_string))
    if(tmp.isEmpty) "none" else tmp.get.substring(tmp.get.indexOf(gene_id_string) + 3)
  }

  /**
    *
    * @param chrm
    * @param sample
    * @param id
    * @param start
    * @param end
    */
  case class ORFseq(chrm: String, sample: String, description: String, id: Int, start: Int, end: Int)

  /**
    * Case class for storing an intergenic sequence as an object
    *
    * @param chrm
    * @param sample
    * @param id
    * @param size
    * @param coords
    */
  case class IntergenicSeq(chrm: String, sample: String, description: String, id: Int, size: Int,
                           coords: Seq[(Int, Int)])


  /**
    * Curried-method to parse a GFFLine collection and return a 4-tuple of a sequence of all (merged)
    * ORFs, the global id for next call, all (valid) intergenic sequences, and the global id for next call.
    *
    *                                Curried parameters:
    * @param isCircular              Molecule is circular
    * @param overlap_proportion      Proportion to use to compute maximum tolerable overlap of ORFs
    * @param minimum_intergenic_size Minimum size required for an intergenic sequence
    * @param showWarnings            Log-level information
    *                                Regular parameters:
    * @param fasta_entry             Fasta-entry of a contig sequence
    * @param sample_name             Name of the sample
    * @param global_orf_id           Starting ORF id to start assigning
    * @param global_inter_id         Starting intergenic id to start assigning
    * @param annotations             Sorted-collection of GFFLine
    * @return (Seq[ORFseq], Int, Seq[(IntergenicSeq)], Int)
    */
  def parseORFs(isCircular: Boolean, overlap_proportion: Double, minimum_intergenic_size: Int, showWarnings: Boolean)
               (fasta_entry: atk.Record,
                sample_name: String,
                global_orf_id: Int,
                global_inter_id: Int,
                annotations: Seq[GFFLine]): (Seq[ORFseq], Int, Seq[(IntergenicSeq)], Int) = {
    //get the name of the contig
    val chrm_name = fasta_entry.getDescription.substring(1).split("\\s+").head
    //get the size of the contig
    val contig_size = fasta_entry.getSequence.size

    /**
      * Create a description for an ORF or an intergenic sequence in the format of:
      * <chrm name>_<sample name>_<some given name>_<start>_<size>
      *
      * @param name  Given name (i.e. gene name)
      * @param start Start coordinate
      * @param end   End coordinate
      * @return String
      */
    def createDescription(name: Option[String], start: Int, end: Int): String = {
      Seq(sample_name, chrm_name, if (name.isEmpty) "" else name.get, start, ((end - start) + 1)).mkString("_")
    }


    /**
      * Method to create an object for an intergenic region. Will check if the intergenic region
      * passes minimum size threshold. Returns a possible intergenic sequence.
      *
      * @param start_position     Last added ORF in the molecule
      * @param first_orf_coords    Ending position of very first ORF
      * @param end_position Position where the ORF should end
      * @param id           Global ID to give
      * @return Option(IntergenicSeq)
      */
    def createIntergeniSeq(start_position: Int, end_position: Int,
                           first_orf_coords: (Int,Int), id: Int): Option[IntergenicSeq] = {
      //get coordinates of potential intergenic sequence aware that the molecule can be circular
      val coordinates = {
        val tmp = Seq((start_position, end_position))
        //for when molecule is circular and reached "end" of the annotations
        if (isCircular && end_position == contig_size) tmp.:+((1, first_orf_coords._1 - 1)) else tmp
      }
      //compute size of intergenic sequence
      val intergenic_size = coordinates.map(x => (x._2 - x._1) + 1).sum
      //check if intergenic region meets minimium size threshold
      if (intergenic_size < minimum_intergenic_size) None
      //meets minimum threshold, create intergenic
      else {
        //create description for intergenic sequence
        val description = createDescription(Option("intergenic"), coordinates.head._1, coordinates.head._1 +
          intergenic_size - 1)
        Option(new IntergenicSeq(chrm_name, sample_name, description, id, intergenic_size, coordinates))
      }
    }

    /** Method to determine whether two GFF objects overlap */
    def isOverlap(x: GFFLine, y: GFFLine): Boolean = {
      //compute the size of the overlap size threshold
      val overlap_threshold = List((x.end - x.start) + 1, (y.end - y.start) + 1).min * overlap_proportion
      //overlap of the two orfs reaches the threshold
      val overlap = min(x.end, y.end) - max(x.start, y.start)
      //check if there is indeed an overlap size of at least the computed threshold
      overlap >= 0 && overlap >= overlap_threshold
    }

    /**
      * Function to determine features of an ORF. If it's a single orf, use as is. If not, there are overlapping ORFs
      * so artificially create an ORF with the max boundaries (smallest start, largest end)
      *
      * @return ORFseq
      */
    def createORFseq(orfs: List[GFFLine], id: Int): ORFseq = {
      //check if there are multiple orfs (overlapping orfs)
      orfs.size match {
        //no overlapping orfs
        case 1 => {
          //create description
          val description = createDescription(orfs.head.name, orfs.head.start, orfs.head.end)
          new ORFseq(chrm_name, sample_name, description, id, orfs.head.start, orfs.head.end)
        }
        //overlapping orfs, create maximal orf
        case _ => {
          //get the largest orf
          val largest_orf = orfs.maxBy(x => (x.end - x.start) + 1)
          //get left-most coordinate
          val min_start = orfs.map(_.start).min
          //get right-most coordinate
          val max_end = orfs.map(_.end).max
          //create description based on largest
          val description = createDescription(largest_orf.name, largest_orf.start, orfs.head.end)
          //log message
          if (showWarnings) println(timeStamp + "----WARNING: overlapping ORFs detected. Using the following maximal start," +
            "end positions: " + (min_start, max_end))
          new ORFseq(chrm_name, sample_name, description, id, largest_orf.start, largest_orf.end)
        }
      }
    }

    /**
      * Tail-recusrive method to iterate through a sequence of gff lines and output the same 4-tuple as the parent
      * method.
      *
      * @param orf_id                Current ORF ID to assign
      * @param inter_id              Current intergenic ID to assign
      * @param remaining_annotations Remaining annotations to go through
      * @param overlaps              The previous annotation along with it's overlapping annotations
      * @param orfs                  Collection with all the (merged) ORFs
      * @param intergenics           Collection with all the (valid) intergenic sequence
      * @return (Seq[ORFseq], Int, Seq[(IntergenicSeq)], Int)
      */
    @tailrec def _parseORFs(orf_id: Int,
                            inter_id: Int,
                            remaining_annotations: Seq[GFFLine],
                            overlaps: List[GFFLine],
                            orfs: Seq[ORFseq],
                            intergenics: Seq[IntergenicSeq]): (Seq[ORFseq], Int, Seq[(IntergenicSeq)], Int) = {
      //no more annotations to go through
      if (remaining_annotations.isEmpty) {
        //there are also no overlapping gene annotations
        if (overlaps.isEmpty) {
          //attempt to create an intergenic sequence
          val intergenic = createIntergeniSeq(orfs.last.end + 1, contig_size, (orfs.head.start, orfs.head.end), inter_id)
          //could not create intergenic sequence, return collections
          if (intergenic.isEmpty) (orfs, orf_id, intergenics, inter_id)
          //return remaining intergenic sequence
          else (orfs, orf_id, intergenics.:+(intergenic.get), inter_id + 1)
        }
        //there are overlapping gene annotations to go through
        else {
          //create a merged orf
          val merged_orfseq = createORFseq(overlaps, orf_id)
          //attempt to create intergenic sequence
          val intergenic = createIntergeniSeq(merged_orfseq.end + 1, contig_size, (orfs.head.start, orfs.head.end), inter_id)
          //add to collection
          if (intergenic.isEmpty) (orfs.:+(merged_orfseq), orf_id + 1, intergenics, inter_id)
          else (orfs.:+(merged_orfseq), orf_id + 1, intergenics.:+(intergenic.get), inter_id + 1)
        }
      }
      //there are still annotations to go through
      else {
        //get next annotation
        val current_annotation = remaining_annotations.head
        //add to overlaps if it is empty or if there is indeed an overlap with any of the previous annotations
        if (overlaps.isEmpty || overlaps.exists(isOverlap(_, current_annotation)))
          _parseORFs(orf_id, inter_id, remaining_annotations.tail, overlaps.:+(current_annotation), orfs, intergenics)
        //no overlap, create orf sequence
        else {
          val orf_seq = createORFseq(overlaps, orf_id)
          //attempt to create intergenic sequence
          val intergenic = {
            //currently in the very first orf
            if(orfs.isEmpty) {
              //genome is circular, wait to create intergenic sequence until the very end
              if(isCircular) None
              //in a genome that is linear, create the intergenic sequence
              else createIntergeniSeq(1, orf_seq.start - 1, (orf_seq.start, orf_seq.end), inter_id)
            }
            //create intergenic seq that is somewhere in the genome
            else createIntergeniSeq(orfs.last.end + 1, orf_seq.start - 1, (orfs.head.start, orfs.head.end), inter_id)
          }
          //add only orf to collection
          if (intergenic.isEmpty) _parseORFs(orf_id + 1, inter_id, remaining_annotations.tail,
            List(current_annotation), orfs.:+(orf_seq), intergenics)
          //add both orf and intergenic
          else _parseORFs(orf_id + 1, inter_id + 1, remaining_annotations.tail,
            List(current_annotation), orfs.:+(orf_seq), intergenics.:+(intergenic.get))
        }
      }
    }
    //fetch all orfs and intergenic sequences
    _parseORFs(global_orf_id, global_inter_id, annotations, List(), Seq(), Seq())
  }

}