package utilities

import scala.util.hashing.MurmurHash3

/**
  * Author: Alex N. Salazar
  * Created on 2-7-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait HashingUtils {

  val zero_byte = 0.toByte
  val one_byte = 1.toByte
  val two_byte = 2.toByte
  val three_byte = 3.toByte
  val min1_byte = (-1).toByte

  /**
    * Method to map nucleotides to bytes.
    *
    * @param c Nucleotide
    * @return Byte
    */
  def encode(c: Char): Byte = c match {
    case 'a' | 'A' => zero_byte
    case 'c' | 'C' => one_byte
    case 'g' | 'G' => two_byte
    case 't' | 'T' => three_byte
    case 'n' | 'N' => min1_byte
    case _ => {
      assert(false, "Cannot encode " + c);
      min1_byte
    }
  }

  def decode(c: Byte): Char = c match {
    case `zero_byte` => 'A'
    case `one_byte` => 'C'
    case `two_byte` => 'G'
    case `three_byte` => 'T'
    case `min1_byte` => 'N'
    case _ => {
      assert(false, "Cannot decode " + c)
      '?'
    }
  }

  /**
    * Method to obtain complement of a nt encoded as byte array
    * @param c Nucleotide encoded as a byte array
    * @return Byte
    */
  def getComplement(c: Byte): Byte =  c match {
    case `zero_byte` => three_byte
    case `one_byte` => two_byte
    case `two_byte` => one_byte
    case `three_byte` => zero_byte
    case `min1_byte` => min1_byte
    case _ => {
      assert(false, "Cannot reverse complement " + c);
      min1_byte
    }
  }

  /**
    * Method to hash byte array using murmurhash
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


}
