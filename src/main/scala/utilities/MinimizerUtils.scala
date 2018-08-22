package utilities

import utilities.FileHandling.timeStamp
import utilities.SequenceUtils.getComplement

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
      val reverse = (getKmerByteHash(kmer.reverse.map(getComplement(_))), 1)
      //ignore if strands are uambiguous
      if (forward._1 == reverse._1) None
      else if (forward._1 < reverse._1) Option(forward)
      else Option(reverse)
    }

    /**
      * Get minimizers given a sequence and a window. Return sequence of hashed kmer value and starting position. The
      * starting position is based on the forward strand even if the kmer originates from the reverse strand.
      *
      * @return Seq[(Int,Int)]
      */
    def getMinimizers(kmer_size: Int, window_size: Int, sequence: Array[Byte]): Set[Minimizer] = {
      val seq_size = sequence.size
      //create an iterator of windows each containing all kmers. w_index is the index of the current window
      sequence.sliding(kmer_size).sliding(window_size).foldLeft((Set[Minimizer](), 1)) {
        case ((minimizers, w_index), window) => {
          //iterate through each kmer in the current window, get hashvalue as well as starting position
          val local_minimizers = window.foldLeft((Seq[Minimizer](), 0)) { case ((local_minimizers, k_index), kmer) => {
            //attempt to get minimum kmer
            val minkmer = getMinKmer(kmer)
            (if (minkmer == None) local_minimizers
            else {
              //compute 1-based position of kmer on forward strand
              val position = {
                if (minkmer.get._2 == 0) w_index + k_index
                else reverseComplementCoord(w_index + k_index, seq_size, kmer_size)
              }
              //create minimizer object
              local_minimizers.:+(new Minimizer(minkmer.get._1, position, minkmer.get._2))
            }
              , k_index + 1)
          }
            //get lowest kmer based on hashvalue and position
          }._1
          if (local_minimizers.isEmpty) {
            println(timeStamp + "--Warning: no minimimizers for read of length: " + sequence.size)
            (minimizers, w_index)
          } else {
            //update minimizer and increment window index
            (minimizers + (local_minimizers.minBy(x => (x.hashvalue, x.orientation))), w_index + 1)
          }
        }
      }._1
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
  }

}
