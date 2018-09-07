package utilities

import java.io.File

import atk.FastaIterator
import utilities.FileHandling.timeStamp
import utilities.SequenceUtils.{getByteComplement, encodeSequence}

import scala.annotation.tailrec
import scala.util.hashing.MurmurHash3

/**
  * Author: Alex N. Salazar
  * Created on 2-7-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object MinimizerUtils {

  /**
    * Minimizer case class
    *
    * @param hashvalue
    * @param position
    * @param orientation 0=forward, 1=reverse
    */
  case class Minimizer(hashvalue: Int, position: Int, orientation: Int)

  trait MMethods {

    /**
      * Function to determine orientation of two minimizer hit orientations:
      * 0 -> node and read are in the same orientation
      * 1 -> node and read are in different orientations
      *
      * @return
      */
    def getOrientation: (Int, Int) => Int = (x, y) => if (x != y) 1 else 0


    /**
      * Generic method to compute the forward-strand position given position, length and kmer size
      *
      * @param i Position
      * @param l Length
      * @param k Kmer size
      * @return Int
      */
    def reverseComplementCoord(i: Int, l: Int, k: Int): Int = l - ((k - 2) + i)

    /**
      * Method to obtain the forward-strand position of a minimizer given the length of the sequence and kmer size
      *
      * @param m Minimizer
      * @param l Length of sequence
      * @param k Kmer size
      * @return Int
      */
    def forwardPosition(m: Minimizer, l: Int, k: Int): Int = {
      if (m.orientation == 0) m.position else reverseComplementCoord(m.position, l, k)
    }

    /**
      * Method to hash byte array using murmurhash
      *
      * @return Long
      */
    def getKmerByteHash: Array[Byte] => Int = kmer => MurmurHash3.bytesHash(kmer)

    /**
      * Get minimum of a kmer. This define as the minimum value of the hash value of a given kmer and it's reverse
      * compliment. Returns tuple of hash value and orienation (0:forward, 1:reverse complement).
      *
      * @return (Int,Int)
      */
    def getMinKmer: Array[Byte] => Option[(Int, Int)] = kmer => {
      //hash forward
      val forward = (getKmerByteHash(kmer), 0)
      //hash reverse
      val reverse = (getKmerByteHash(kmer.reverse.map(getByteComplement(_))), 1)
      //ignore if strands are uambiguous
      if (forward._1 == reverse._1) None
      else if (forward._1 < reverse._1) Option(forward)
      else Option(reverse)
    }

    /**
      * Get minimizers given a sequence and a window. Return sequence of hashed kmer value and starting position.
      *
      * @param kmer_size      Size of kmers
      * @param window_size    Window size for minimizer
      * @param sequences      Byte-encoded sequence (output of SequenceUtils.encodeSequence function)
      * @param total_seq_size Total size of sequence
      * @return Set[Minimizer]
      */
    def getMinimizers(kmer_size: Int,
                      window_size: Int,
                      sequences: Array[(Array[Byte], Int)],
                      total_seq_size: Int): List[Minimizer] = {
      //iterate through each encoded sequence and collect minimizers
      sequences.foldLeft(List[Minimizer]()){ case (global_minimizers, (sequence, starting_position)) => {
        //get sequence size
        val seq_size = sequence.size
        //check that that the sequence is at least the same size of desired kmers
        if (seq_size < kmer_size) global_minimizers
        else {
          //create an iterator of windows each containing all kmers. w_index is the index of the current window
          sequence.sliding(kmer_size).sliding(window_size).foldLeft((global_minimizers, 0)) {
            case ((minimizers, w_index), window) => {
              //iterate through each kmer in the current window, get hashvalue as well as starting position
              val local_minimizers = window.foldLeft((List[Minimizer](), 0)) {
                case ((local_minimizers, k_index), kmer) => {
                  //attempt to get minimum kmer
                  val minkmer = getMinKmer(kmer)
                  (if (minkmer == None) local_minimizers
                  else {
                    //compute 1-based position of kmer corresponding strand
                    val position = {
                      if (minkmer.get._2 == 0) w_index + k_index + starting_position
                      else reverseComplementCoord(w_index + k_index + starting_position, total_seq_size, kmer_size)
                    }
                    //create minimizer object
                    (new Minimizer(minkmer.get._1, position, minkmer.get._2)) :: local_minimizers
                  }
                    , k_index + 1)
                }
              }._1
              if (local_minimizers.isEmpty) {
                println(timeStamp + "--Warning: no minimimizers for read of length: " + sequence.size)
                (minimizers, w_index)
              } else {
                //update minimizer and increment window index
                ((local_minimizers.minBy(x => (x.hashvalue, x.orientation))) :: minimizers, w_index + 1)
              }
            }
          }._1
        }
      } //remove duplicates
      }.groupBy(x => (x.hashvalue, x.orientation, x.position)).toList.map(_._2.head)
    }

    /**
      * Method to cluster minimizer hits. Returns clusters of minimizers hits as 2-tuple (minimizer, forward position).
      *
      * @param minimizers
      * @param seq_size Size of read
      * @param max_dist
      * @return Seq[Set[(MinKmer,Int)]
      **/
    def clusterMinimizers(minimizers: List[Minimizer], kmer_size: Int, seq_size: Int,
                          max_dist: Int, debug: Boolean): Seq[Set[(Minimizer, Int)]] = {

      /**
        * Tail-recursive implementation of clustering procedure
        *
        * @param remaining
        * @param last Last position observed in forward strand
        * @param acc
        * @param clusters
        * @return
        */
      @tailrec def _clusterMinimizers(remaining: List[(Minimizer, Int)], last: Int, acc: Set[(Minimizer, Int)],
                                      clusters: Seq[Set[(Minimizer, Int)]]): Seq[Set[(Minimizer, Int)]] = {
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
          if (x.orientation == 0) x.position else reverseComplementCoord(x.position, seq_size, kmer_size)
        }
        (x, forward_position)
      }).sortBy(_._2)
      if (debug) println(timeStamp + "----Clustering the following minimizers: " + forward_map_minimizers)
      _clusterMinimizers(forward_map_minimizers, minimizers.head.position, Set(), Seq())
    }

    /**
      * Double curried-method to fetch minimizers from either ORF or Intergenic sequences. First set of parameters
      * are just standard configurations. Second is an Option[Map] which determines whether minimizers will be
      * extracted from orf (provide a Map[orf id -> node id] or intergenic sequences (empty map).
      *
      * For ORF minimizers, creates and returns a minimizer hashtable and node info hashtable:
      * The minimizer hashtable schema: kmer -> Map(node ID -> Set(position, orientation))
      * The node info hashtable schema: node id -> Map(orf id -> (minimizer set size, seq size))
      *
      * For Intergenic minimizers, same as above but the minimizer table only contains keys and
      * the values are empty while the node info is empty.
      *
      * In both cases, returns minimizer hashtable and node info hashtable.
      *
      * Configuration parameters:
      *
      * @param kmer_size    Size of kmers
      * @param window_size  Window size for minimizers
      * @param debug        Show debug-level information
      *                     Providing a non-empty maps turns on the ORF minimizer procedure:
      * @param node_map     Map of orf ID -> node ID
      * @param nodeinfo_map Accumulating node info hashtable (see description above)
      *                     Sequence-based parameters:
      * @param fasta_file   Fasta file of some sample
      * @param kmer_index   Accumulating kmer index hashtable (see description above)
      * @return (Map[Int, Map[Int, Set[(Int, Int)], Map[Int, Map[Int, (Int, Int)])
      */
    def fetchMinimizers(kmer_size: Int, window_size: Int, debug: Boolean)
                       (node_map: Option[Map[Int, Int]], nodeinfo_map: Option[Map[Int, Map[Int, (Int, Int)]]])
                       (fasta_file: File, kmer_index: Map[Int, Map[Int, Set[(Int, Int)]]]
                       ): (Map[Int, Map[Int, Set[(Int, Int)]]], Map[Int, Map[Int, (Int, Int)]]) = {
      //create fasta iterator
      val fasta_iterator = new FastaIterator(fasta_file)

      /**
        * Tail-recusive method to iterate through fasta entries and update hashtable for orf sequences
        *
        * @param _kmer_index Hashtable to update
        * @return Map[Int, Map[Int,Int]
        **/
      @tailrec def fetchORFMinimizers(_kmer_index: Map[Int, Map[Int, Set[(Int, Int)]]],
                                      _nodeinfo_map: Map[Int, Map[Int, (Int, Int)]]
                                     ): (Map[Int, Map[Int, Set[(Int, Int)]]], Map[Int, Map[Int, (Int, Int)]]) = {
        if (!fasta_iterator.hasNext) (_kmer_index, _nodeinfo_map)
        else {
          //get current orf
          val current_orf = fasta_iterator.next()
          //get orf id and sequence
          val (orf_id, orf_seq) = (current_orf.getDescription.substring(1).toInt, current_orf.getSequence)
          //get corresponding node id, for intergenic return none
          val node_id = node_map.get(orf_id)
          //get minimizers
          val minimizers = getMinimizers(kmer_size, window_size, encodeSequence(orf_seq), orf_seq.size)
          //get UNIQUE minimizers; note: a hashvalue can belong to multiple minimizers
          val unique_minimizers_size = minimizers.map(_.hashvalue).size
          //update node info map
          val updated_nodeinfo = {
            val current = _nodeinfo_map.getOrElse(node_id, Map[Int, (Int, Int)]())
            assert(current.get(orf_id) == None, "Multiple occurrences of ORF sequence: " + orf_id)
            _nodeinfo_map + (node_id -> (current + (orf_id -> (unique_minimizers_size, orf_seq.size))))
          }
          //update kmer index
          fetchORFMinimizers(minimizers.foldLeft(_kmer_index)((local_kmer_index, minimizer) => {
            //get kmer's map of node ids to starting positions and update accordingly
            val current = {
              //get all nodes where this kmer appears in
              val node_map = local_kmer_index.getOrElse(minimizer.hashvalue, Map[Int, Set[(Int, Int)]]())
              //get all starting positions of current node
              val node_seq = node_map.getOrElse(node_id, Set[(Int, Int)]())
              //update starting positions of current node
              node_map + (node_id -> (node_seq + ((minimizer.position, minimizer.orientation))))
            }
            //update kmer index
            local_kmer_index + (minimizer.hashvalue -> current)
          }), updated_nodeinfo)
        }
      }

      /**
        * Tail-recursive method as above but to get intergenic sequence minimizers
        *
        * @param _kmer_index
        * @return (Map[Int, Map[Int, Set[(Int, Int)], Map[Int, Map[Int, (Int, Int)])
        */
      @tailrec def fetchIntergenicMinimizers(_kmer_index: Map[Int, Map[Int, Set[(Int, Int)]]]
                                            ): (Map[Int, Map[Int, Set[(Int, Int)]]], Map[Int, Map[Int, (Int, Int)]]) = {
        if (!fasta_iterator.hasNext) (_kmer_index, Map())
        else {
          //get current orf
          val current_orf = fasta_iterator.next()
          //get orf id and sequence
          val (inter_id, inter_seq) = (current_orf.getDescription.substring(1).toInt, current_orf.getSequence)
          //get corresponding node id, for intergenic return none
          //val node_id = node_map.get(inter_id)
          //update UNIQUE minimizers; note: a hashvalue can belong to multiple minimizers
          fetchIntergenicMinimizers(
            getMinimizers(kmer_size, window_size, encodeSequence(inter_seq), inter_seq.size)
              .foldLeft(_kmer_index)((b, a) => b + (a.hashvalue -> Map())))
        }
      }
      //call appropriate method when dealing with orfs or intergenic sequences
      if (node_map.isEmpty) fetchIntergenicMinimizers(kmer_index) else fetchORFMinimizers(kmer_index, nodeinfo_map.get)
    }
  }

}
