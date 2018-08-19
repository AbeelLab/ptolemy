package utilities

import utilities.FileHandling.timeStamp

/**
  * Author: Alex N. Salazar
  * Created on 15-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait SequenceUtils extends HashingUtils{

  case class MinKmer(hashvalue: Int, position: Int, orientation: Int)

  /**
    * Get minimizers given a sequence and a window. Return sequence of hashed kmer value and starting position. The
    * starting position is based on the forward strand even if the kmer originates from the reverse strand.
    *
    * @return Seq[(Int,Int)]
    */
  def getMinimizers(kmer_size: Int, window_size: Int, sequence: Array[Byte]): Set[MinKmer] = {
    val seq_size = sequence.size
    //create an iterator of windows each containing all kmers. w_index is the index of the current window
    sequence.sliding(kmer_size).sliding(window_size).foldLeft((Set[MinKmer](), 1)) {
      case ((minimizers, w_index), window) => {
      //iterate through each kmer in the current window, get hashvalue as well as starting position
      val local_minimizers = window.foldLeft((Seq[MinKmer](), 0)) { case ((local_minimizers, k_index), kmer) => {
        //attempt to get minimum kmer
        val minkmer = getMinKmer(kmer)
        (if (minkmer == None) local_minimizers
        else {
          //compute 1-based position of kmer on forward strand
          val position = {
            if(minkmer.get._2 == 0) w_index + k_index
            else reverseComplementCoord(w_index + k_index, seq_size, kmer_size)
          }
          //create minimizer object
          local_minimizers.:+(new MinKmer(minkmer.get._1, position, minkmer.get._2))
        }
        , k_index + 1)
      }
        //get lowest kmer based on hashvalue and position
      }._1
        if(local_minimizers.isEmpty){
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
    * Method to compute the corresponding reverse complement coordinate for a given coordinate.
    * @param i Position
    * @param l Length of sequence
    * @param k Kmer size
    * @return Int
    */
  def reverseComplementCoord(i: Int, l: Int, k: Int): Int = l - ((k - 2) + i)

}
