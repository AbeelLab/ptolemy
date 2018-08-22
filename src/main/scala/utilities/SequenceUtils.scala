package utilities

import utilities.FileHandling.timeStamp

import scala.annotation.tailrec
import scala.util.hashing.MurmurHash3

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

  /**
    * Method to decode an array of bytes to DNA sequence
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
      assert(false, "Cannot decode " + c)
      '?'
    }
  }

  /**
    * Method to obtain complement of a nt encoded as byte array
    *
    * @param c Nucleotide encoded as a byte array
    * @return Byte
    */
  def getComplement(c: Byte): Byte = c match {
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


  trait FastaUtils {

    /**
      * Case class for FASTA-sequence entry. Sequence is represented as a byte array
      *
      * @param name
      * @param sequence
      */
    case class FastaEntry(name: String, sequence: Array[Byte])


    /**
      * Method to load set a chunk of reads up to some memory usage.
      *
      * @param iterator
      * @param acc_iterator
      * @param acc_name
      * @param acc_sequence
      * @return List[FastaEntry]
      */
    def loadFastaChunk(iterator: Iterator[String],
                       acc_iterator: List[FastaEntry],
                       acc_nreads: Int,
                       minsize: Int,
                       acc_name: Option[String],
                       acc_sequence: StringBuilder): (List[FastaEntry], Option[String]) = {
      //empty iterator
      if (!iterator.hasNext) {
        //no additional reads means reached end of file, check that there is a remaining read with sequence
        assert(acc_name != None, "Expected FASTA-entry at end of file")
        assert(!acc_sequence.isEmpty, "Expected corresponding sequence at end of file")
        //add last read to iterator
        (acc_iterator.:+(new FastaEntry(acc_name.get.substring(1), acc_sequence.toArray.map(encode(_)))), None)
      }
      //loaded max number of reads
      else if (acc_nreads == 0) {
        //check that there is a read entry in the acc
        assert(!acc_name.isEmpty, "Expected read entry after loading max number of reads: " + iterator.next())
        //check that the corresponding sequence is empty
        assert(acc_sequence.isEmpty, "Expected no corresponding sequence for read entry after loading max number of " +
          "reads: " + acc_sequence.mkString(""))
        (acc_iterator, acc_name)
      }
      //keep loading reads
      else {
        //get the next line
        val line = iterator.next()
        //for when the FASTA acc entry is empty, line is expected to be a new FASTA entry (e.g. beginning of iteration)
        if (acc_name == None) {
          //sanity checks
          assert(line.startsWith(">"), "Expected start of a FASTA sequence: " + line)
          assert(acc_sequence.isEmpty, "Expected empty sequence accumulator for " + line)
          //add as start of new FASTA entry
          loadFastaChunk(iterator, acc_iterator, acc_nreads, minsize, Option(line), acc_sequence)
        }
        else {
          //sanity check that the next line starts with entry character and sequence is empty
          assert(acc_name.get.startsWith(">"), "Expected line to be start of a FASTA-entry: " + acc_name.get)
          //for when there is an existing FASTA-entry and the next line is a new FASTA-entry
          if (line.startsWith(">")) {
            //sanity check
            assert(acc_sequence != None, "Expected FASTA sequence for entry " + acc_name.get + ", instead" +
              " found start of new entry: " + line)
            //only add if fasta entry meets minimum read length
            if (acc_sequence.length < minsize)
              loadFastaChunk(iterator, acc_iterator, acc_nreads, minsize, Option(line), new StringBuilder)
            else {
              //create a new FastaEntry object
              val fasta_entry = new FastaEntry(acc_name.get.substring(1), acc_sequence.toArray.map(encode(_)))
              //add accordingly
              loadFastaChunk(iterator, acc_iterator.:+(fasta_entry), acc_nreads - 1, minsize, Option(line), new StringBuilder)
            }
          }
          //for when there is an existing FASTA entry and the next line is a sequence line
          else loadFastaChunk(iterator, acc_iterator, acc_nreads, minsize, acc_name, acc_sequence.append(line))
        }
      }
    }
  }

}
