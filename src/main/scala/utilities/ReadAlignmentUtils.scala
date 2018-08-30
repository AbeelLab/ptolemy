package utilities

import java.io.File

import utilities.FileHandling.{openFileWithIterator, timeStamp}
import utilities.MinimizerUtils.{MMethods, Minimizer}
import utilities.NumericalUtils.{abs, max, min}
import utilities.GFAutils.ReadGFA
import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 25-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait ReadAlignmentUtils extends MMethods with ReadGFA with ReadUtils {

  /**
    * Curried method to fetch minimizers from a read. First set of parameter are configuration parameters. Second set
    * is the read as a FastaEntry object. Returns minimizers as a map of kmer -> sequence(minimizers)
    *
    * @param kmer_size   Size of kmer
    * @param window_size Window size
    * @param read        FastaEntry object of a read
    * @return Map[Int, Seq[MinKmer]
    */
  def getReadMinimizers(kmer_size: Int, window_size: Int)
                       (read: FastaEntry): (Map[Int, Seq[Minimizer]]) = {
    //first get minimizers
    getMinimizers(kmer_size, window_size, read.sequence, read.length)
      //then build map, note that the same hashvalue can have multiple minimizers
      .foldLeft(Map[Int, Seq[Minimizer]]())((map, minimizer) => {
      //get current minimizers for current hashvalue
      val current = map.getOrElse(minimizer.hashvalue, Seq[Minimizer]())
      map + (minimizer.hashvalue -> (current.:+(minimizer)))
    })
  }


  /**
    * Curried-method for predicting the start and end coordinates of a node on a read. First set of parameters are
    * configuration parameters per read. Second set of parameters are per node.
    *
    * @param read_name    Read ID
    * @param debug        Log-level
    * @param node         Node ID
    * @param orientation  Orientation of the node
    * @param avg_seq_size Size of the node
    * @param cluster      Sequence of 2-tuples corresponding to the minimizer hits as (node position, read position).
    *                     Assumes ascending order on read position.
    * @return (Int,Int)
    */
  def predictCoords(read_name: String, debug: Boolean)
                   (node: Int, orientation: Int, avg_seq_size: Int, cluster: Seq[(Int, Int)]): (Int, Int) = {
    /**
      * Tail-recursive method to go through the cluster and find the most logical start/end coordinates.
      *
      * @param sorted
      * @return (Int,Int)
      */
    @tailrec def _predictCoords(sorted: Seq[(Int, Int)], left: Boolean): (Int, Int) = {
      if (sorted.size <= 1) (-1, -1)
      else {
        //get left most position
        val left_most = sorted.head
        //get right most position
        val right_most = sorted.last
        //compute predicted star/end coordinates on the read
        val (start, end) = {
          //calculate according to the orientation
          if (orientation == 0) (left_most._1 - left_most._2, right_most._1 + abs(avg_seq_size - right_most._2))
          else (right_most._1 - abs(avg_seq_size - right_most._2), left_most._1 + left_most._2)
        }
        //found a logical start/end coordinate
        if (start < end) (start, end)
        //try again using the next left/right-most positions
        else {
          if (debug) println(timeStamp + "----WARNING: predicted start/end is illogical: " + (left_most, right_most))
          //determine direction to subset
          // take away from left
          if (left) _predictCoords(sorted.tail, false)
          // take away form right
          else _predictCoords(sorted.take(sorted.size - 1), true)
        }
      }
    }

    _predictCoords(cluster, true)
  }

  /**
    * Method to curate a given read alignment by validating aigned nodes and extracting unaligned regions in the read
    *
    * @param alignment          Sorted sequence of aligned nodes as a 4-tuple (node, (start,end), orientation, hits)
    * @param read_size          Size of read
    * @param min_size           Minimum unaligned region size threshold
    * @param overlap_proportion Proportion of smallest node alignment to be used as the maximum tolerable overlap
    *                           between two nodes
    * @return
    */
  def curateAlignments(alignment: Seq[(Int, (Int, Int), Int, Int)],
                       read_size: Int,
                       min_size: Int,
                       overlap_proportion: Double,
                      ): Seq[(Int, (Int, Int), Int, Int)] = {

    /**
      * Method to determine whether two node alignments overlap by some threshold. If so, these are considered invalid
      * alignments.
      *
      * @param x Start and end coordinates of a node alignment
      * @param y Same as x
      * @return Boolean
      */
    def isInvalidAlignment(x: (Int, Int), y: (Int, Int)): Boolean = {
      //compute the size of the overlap size threshold
      val overlap_threshold = (min((x._2 - x._1) + 1, (y._2 - y._1) + 1) * overlap_proportion).toInt
      //overlap of the two orfs reaches the threshold
      val overlap = min(x._2, y._2) - max(x._1, y._1)
      //check if there is indeed an overlap size of at least the computed threshold
      overlap > overlap_threshold
    }

    /**
      * Method to merge invalid alignments into an ambiugous representation
      *
      * @param m Overlapping invalid alignments
      * @return (Int, (Int, Int), Int, Int)
      */
    def mergeAlignments(m: Seq[(Int, (Int, Int), Int, Int)]): (Int, (Int, Int), Int, Int) = {
      //get left-most coordinate
      val start = m.map(_._2._1).min
      //get right-most coordinate
      val end = m.map(_._2._2).max
      //return ambiguous representation
      (-1, (start, end), -1, -1)
    }

    /**
      * Tail-recursive method to merge invalid node alignments into a single, ambiguous representation
      *
      * @param a      Initial sequence of alignments in ascending order
      * @param acc    Accumulating collection of invalid alignments
      * @param merged Final sequence of alignments with the ambiguous representations
      * @return Seq[(Int, (Int, Int), Int, Int)]
      */
    @tailrec def mergeInvalids(a: Seq[(Int, (Int, Int), Int, Int)],
                               acc: Seq[(Int, (Int, Int), Int, Int)],
                               merged: Seq[(Int, (Int, Int), Int, Int)]): Seq[(Int, (Int, Int), Int, Int)] = {
      //no more nodes to process
      if (a.isEmpty) {
        //no remaining nodes to merge
        if (acc.isEmpty) merged
        //merge remaining nodes if there is more than 1
        else if (acc.size == 1) merged.:+(acc.head) else merged.:+(mergeAlignments(acc))
      }
      //accumulative is empty so add current node alignment
      else if (acc.isEmpty) mergeInvalids(a.tail, acc.:+(a.head), merged)
      //check if current node alignment overlaps/is invalid with any current invalid alignments
      else {
        //current node alignment overlaps with one of the invalid node alignments
        if (acc.exists(x => isInvalidAlignment(x._2, a.head._2))) mergeInvalids(a.tail, acc.:+(a.head), merged)
        //current node alignment does not overlap, merge current overlapping one and move to collection
        else mergeInvalids(a.tail, Seq(a.head),
          if (acc.size == 1) merged.:+(acc.head) else merged.:+(mergeAlignments(acc)))
      }
    }

    //merge invalid alignments
    val merged_alignments = mergeInvalids(alignment, Seq(), Seq())
    //compute ending index
    val ending_index = if (merged_alignments.size == 1) 0 else (merged_alignments.size - 2 + 1) - 1

    /**
      * Curate merged alignments: create a sequence of unaligned, invalid, and valid alignments.
      */
    merged_alignments.sliding(2).foldLeft((0, Seq[(Int, (Int, Int), Int, Int)]())) {
      case ((index, curated_alignments), nodes) => {
        /**
          * Get unaligned region between the first and second node
          *
          * @return (Int,Int)
          */
        def getUnalignedRegion(): (Int, Int) = (nodes.head._2._2 + 1, nodes(1)._2._1 - 1)

        /**
          * Create alignment instance of an unaligned region. This is determined by the first element of the 4-tuple
          * which is -2 flag
          *
          * @return
          */
        def unalignedInstance: ((Int, Int)) => (Int, (Int, Int), Int, Int) = x => (-2, (x._1, x._2), -1, -1)

        val local_curation = {
          //construct according to index
          index match {
            //first and last index: add the following
            case x if (x == 0 && x == ending_index) => {
              //this is to handle cases when there is only one node alignment: first add common regions
              val tmp = {
                //beginning
                Seq(unalignedInstance((1, nodes.head._2._1 - 1)),
                  //first alignment
                  nodes(0))
              }
              //if tehre is only one alignment, add end of read
              if (nodes.size == 1) tmp.:+(unalignedInstance((nodes(0)._2._2 + 1, read_size)))
              //multiple alignments, add the following
              else tmp ++
                Seq(//in-between
                  unalignedInstance(getUnalignedRegion),
                  //second alignment
                  nodes(1),
                  //end of read
                  unalignedInstance((nodes(1)._2._2 + 1, read_size)))
            }
            //very first index: add the following
            case 0 =>
              //beginning
              Seq(unalignedInstance((1, nodes.head._2._1 - 1)),
                //first alignment
                nodes(0),
                //in-between
                unalignedInstance(getUnalignedRegion()))
            //last index: add the following
            case `ending_index` =>
              //first alignment
              Seq(nodes(0),
                //in-between
                unalignedInstance(getUnalignedRegion()),
                //last alignment
                nodes(1),
                //end of read
                unalignedInstance((nodes(1)._2._2 + 1, read_size)))
            //any other index, add the following
            case _ =>
              //first alignment and unaligned instance
              Seq(nodes(0), unalignedInstance(getUnalignedRegion()))
          }
        }
        //increment index, add local curations while filter for unaligned regions (flag -2) meeting minimum size
        (index + 1, curated_alignments ++ local_curation.filter(x => x._1 != -2 || x._2._2 - x._2._1 + 1 >= min_size))
      }
    }._2
  }

  /**
    * Function to obtain the sequence of maximum contiguous alignments given a sequence of node alignments. Returns a
    * sequence of maximum contiguous alignments as 2-tuple (node id, orientation).
    *
    * @return Seq[Seq[(Int,Int)]
    */
  def maxContiguousAlignments: Seq[(Int, (Int, Int), Int, Int)] => Seq[Seq[(Int, Int)]] = alignments => {

    /**
      * Tail-recursive method to find maximum contiguous alignments. Returns same as parent function.
      *
      * @param _alignments      Node alignments
      * @param acc              Accumulating maximum contiguous alignment
      * @param split_alignments Sequence of split alignments
      * @return Seq[Seq[(Int, Int)]
      */
    @tailrec def _maxContiguousAlignments(_alignments: Seq[(Int, (Int, Int), Int, Int)],
                                          acc: Seq[(Int, Int)], split_alignments: Seq[Seq[(Int, Int)]]): Seq[Seq[(Int, Int)]] = {
      //not more alignments to go through
      if (_alignments.isEmpty) {
        //no max contiguous alignment to add, or add max contiguous alignment
        if (acc.isEmpty) split_alignments else split_alignments.:+(acc)
      }
      else {
        //get current alignment
        val current = (_alignments.head._1, _alignments.head._3)
        //if it's not an ambiguous alignment or unaligned region, add to current maximum contiguous alignment
        if (current._1 >= 0) _maxContiguousAlignments(_alignments.tail, acc.:+(current), split_alignments)
        //ambiguous alignment, create a split alignment
        else _maxContiguousAlignments(_alignments.tail, Seq(),
          if (acc.isEmpty) split_alignments else split_alignments.:+(acc))
      }
    }

    _maxContiguousAlignments(alignments, Seq(), Seq())
  }


  /**
    * Function to determine the orientation of an alignment. Note: this is simple majority vote.
    * '+' -> Forward
    * '-' -> Reverse
    * '?' -> Ambiguous
    *
    * @return Char
    */
  def determineAlignmentOrientation: Seq[(Int, Int)] => Char = alignments => {
    /**
      * There should only be two orientation options, forward (0) or reverse (1)
      *
      * @return Boolean
      */
    def isValid: Int => Boolean = x => x == 0 || x == 1

    //compute majority vote
    val tmp = alignments.foldLeft((0, 0)) { case ((f, r), alignment) => {
      assert(isValid(alignment._2), "Invalid orientation: " + alignment)
      if (alignment._2 == 0) (f + 1, r) else (f, r + 1)
    }
    }
    // /majority vote, else ambiguous
    if (tmp._1 == tmp._2) '?'
    else if (tmp._1 > tmp._2) '+'
    else '-'
  }

  /**
    * Method to obtain edge and node coverage from a tmp alignment file as well as a few statistics
    *
    * @param file
    * @return
    */
  def getAlignmentStats(file: File): (Map[(Int, Int), Int], Map[Int, Int], Int, Int) = {
    //open file as iterator
    val lines = openFileWithIterator(file)

    /**
      * Tail-recursive method to iterate through tmp alignment file and obtain edge and node coverages as well as
      * some statistics
      *
      * @param edge_coverage
      * @param node_coverage
      * @param total_alignments
      * @param total_unmapped
      * @return
      */
    @tailrec def _getAlignmentStats(edge_coverage: Map[(Int, Int), Int],
                                    node_coverage: Map[Int, Int],
                                    total_alignments: Int,
                                    total_unmapped: Int
                                   ): (Map[(Int, Int), Int], Map[Int, Int], Int, Int) = {
      if (!lines.hasNext) (edge_coverage, node_coverage, total_alignments, total_unmapped)
      else {
        //get next line
        val line = lines.next()
        //don't do anything if it's empty
        if (line.isEmpty) _getAlignmentStats(edge_coverage, node_coverage, total_alignments, total_unmapped)
        else {
          //parse alignment line
          val (read_id, alignment) = parseAlignment(line)
          //update stats
          _getAlignmentStats(
            //if unmapped or maps to only one node, do not update edge coverage
            (if (alignment.size < 2) edge_coverage
            //else updated coverage of observed edges as well as total number of reads
            else alignment.sliding(2).foldLeft(edge_coverage)((local_ecoverage, _edge) => {
              val edge = (_edge(0), _edge(1))
              val current_coverage = local_ecoverage.getOrElse(edge, 0)
              local_ecoverage + (edge -> (current_coverage + 1))
            })),
            //incremenet node coverages
            (alignment.foldLeft(node_coverage)((local_ncoverage, node) => {
              val current_coverage = local_ncoverage.getOrElse(node, 0)
              local_ncoverage + (node -> (current_coverage + 1))
            })),
            //increment total reads and total mapped
            total_alignments + 1, if (alignment.isEmpty) total_unmapped + 1 else total_unmapped
          )
        }
      }
    }

    _getAlignmentStats(Map(), Map(), 0, 0)
  }

}
