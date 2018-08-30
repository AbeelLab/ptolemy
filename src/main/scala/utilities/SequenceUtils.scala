package utilities

import scala.collection.mutable.ArrayBuffer

/**
  * Author: Alex N. Salazar
  * Created on 15-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */


object SequenceUtils {

  /**
    * Encoding bytes
    */
  private val zero_byte = 0.toByte
  private val one_byte = 1.toByte
  private val two_byte = 2.toByte
  private val three_byte = 3.toByte
  private val min1_byte = (-1).toByte

  /**
    * Function to encode a nucleotide into a byte. For standard nucleotides, returns option(byte). Any
    * other returns none
    *
    * @return Option[Byte]
    */
  def encode: Char => Option[Byte] = n => {
    n match {
      case 'a' | 'A' => Option(zero_byte)
      case 'c' | 'C' => Option(one_byte)
      case 'g' | 'G' => Option(two_byte)
      case 't' | 'T' => Option(three_byte)
      case _ => None
    }
  }

  /**
    * Method to encode a given sequence into bytes. It considers IUPAC nucleotides: when encountering a
    * non-standard nucleotide (A|T|C|G), the nucleotide is discarded and the upstream sequence is taken as the
    * longest contiguous encoded sequence. The result is thus an array of arrays representing all longest encoding
    * sequences. Note: it uses a non-functional and mutable approach since it internally uses array buffers.
    *
    * @return Array[Array[Byte]
    */
  def encodeSequence(seq: String): Array[(Array[Byte], Int)] = {
    /**
      * Initiate mutable array buffer of array buffer containing all contiguous encoded sequences along with starting
      * positions
      */
    var byte_arrays = ArrayBuffer[(Array[Byte], Int)]()
    /**
      * Initiate mutable array buffer for accumulating contiguous encoded sequence
      */
    var acc = ArrayBuffer[Byte]()
    /**
      * Initiate starting position for the current accumulating contiguous encodeds sequence
      */
    var starting_position = 1
    /**
      * Initiate starting position for the current position at any point
      */
    var current_position = 1
    //iterate through each nucleotide in the sequence and extract all longest contiguous encoded sequences
    seq.toCharArray.foreach(nt => {
      //attempt to encode current nucleotide
      val encoding = encode(nt)
      //encountered a non-standard nucleotide
      if (encoding.isEmpty) {
        //if acc is not empty, add the upstream sequence to the array of encoded sequences
        if (!acc.isEmpty) {
          byte_arrays += ((acc.toArray, starting_position))
          //empty acc
          acc = ArrayBuffer[Byte]()
        }
        //increment positions regardless
        current_position += 1;
        starting_position = current_position
      }
      //still in a contiguou sequence, add encoded nucleotide
      else {
        //add nucleotide
        acc += encoding.get
        //increment current position
        current_position += 1
      }
    })

    //add remaining contiguous sequence, if it exists
    if (acc.isEmpty) byte_arrays.toArray else (byte_arrays += ((acc.toArray, starting_position))).toArray
  }

  /**
    * Method to decode an array of bytes to DNA sequence
    *
    * @param c
    * @return
    */
  def decode(c: Byte): Char = c match {
    case `zero_byte` => 'A'
    case `one_byte` => 'C'
    case `two_byte` => 'G'
    case `three_byte` => 'T'
    case `min1_byte` => 'N'
    case _ => {
      assert(false, "Cannot decode nucleotide " + c)
      '?'
    }
  }

  /**
    * Function to get the complement for a given nucleotide.
    *
    * @return Byte
    */
  def getByteComplement: Byte => Byte = b => {
    b match {
      case `zero_byte` => three_byte
      case `one_byte` => two_byte
      case `two_byte` => one_byte
      case `three_byte` => zero_byte
      case _ => {
        assert(false, "Cannot obtain complement nucleotide for byte " + b);
        min1_byte
      }
    }
  }
}
