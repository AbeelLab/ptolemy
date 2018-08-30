package utilities

import java.io.File

import utilities.SequenceUtils.{encodeSequence}

import scala.annotation.tailrec
import scala.collection.parallel.ParSeq

/**
  * Author: Alex N. Salazar
  * Created on 28-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

trait ReadUtils {

  val fasta_extension = "FASTA"
  val fastq_extension = "FASTQ"
  val fastq_gz_extension = "FASTQ.gz"

  /**
    * Case class for FASTA-sequence entry. Sequence is represented as a byte array
    *
    * @param name
    * @param sequence
    */
  case class FastaEntry(name: String, sequence: Array[(Array[Byte], Int)], length: Int)

  def determineFileType(file: File): String = {
    /**
      * Function to determine whether given file extension is a FASTA-formatted file
      *
      * @return Option[String]
      */
    def isFasta: Unit => Option[String] = Unit => {
      val t = List("fasta", "fna", "fa").exists(x => file.getName.endsWith(x))
      if (t) Option(fasta_extension) else None
    }

    /**
      * Function to determine whether given file extension is a FASTQ-formatted file
      *
      * @return Option[String]
      */
    def isFastq: Unit => Option[String] = Unit => {
      val t = List("fastq").exists(x => file.getName.endsWith(x))
      if (t) Option(fastq_extension) else None
    }

    /**
      * Function to determine whether given file extension is a FASTQ-GZIPPED-formatted file
      *
      * @return Option[String]
      */
    def isFastqGZ: Unit => Option[String] = Unit => {
      val t = List("fastq.gz").exists(x => file.getName.endsWith(x))
      if (t) Option(fastq_gz_extension) else None
    }

    //empty option
    var file_type: Option[String] = None
    //list of all functions
    val functions = List(isFasta, isFastq, isFastqGZ)
    //iterate through each function and attempt to find file type
    functions.foreach(f => {
      //check whether there is matching file type
      val t = f()
      if (!t.isEmpty) file_type = t
    })
    //sanity check
    assert(!file_type.isEmpty, "Could not infer file type for given reads: " + file.getName)
    //return type
    file_type.get
  }


  /**
    * Function to extract read ID from a given FASTQ or FASTA line
    *
    * @return String
    */
  def getReadName: String => String = line => line.substring(1).split("\\s+").head

  /**
    * Function to convert a list of FASTA or FASTQ lines to FastaEntry object
    *
    * @return FastaEntry
    */
  def fastq2FastaEntry: List[String] => FastaEntry = lines =>
    new FastaEntry(getReadName(lines.head), encodeSequence(lines(1)), lines(1).size)

  /**
    * Tail-recursive method to load some maximum number of reads into a chunk from a FASTQ-formatted file.
    *
    * @param iterator        FASTQ iterator
    * @param acc             Accumulating lines for a single read
    * @param remaining_reads Remaining number of reads to load
    * @param minsize         Minimum read size to load
    * @param chunk           Accumulating chunk of reads
    * @return (List[FastaEntry], Option[String])
    */
  @tailrec final def loadFastqChunk(iterator: Iterator[String],
                                    acc: List[String],
                                    remaining_reads: Int,
                                    minsize: Int,
                                    chunk: List[FastaEntry]): (ParSeq[FastaEntry], Option[String]) = {
    //no more reads to process or loaded maximum number of reads
    if (iterator.isEmpty || remaining_reads == 0) {
      //no accumulating reads left
      if (acc.isEmpty) (chunk.par, None)
      else {
        //sanity check
        assert(acc.size == 4, "Unexpected number of lines in accumulating reads after loading max read chunk: " + acc)
        //great fasta entry
        val read = fastq2FastaEntry(acc)
        //only add remaining read if it meets minimum size threshold
        if (read.length < minsize) (chunk.par, None) else (chunk.:+(fastq2FastaEntry(acc)).par, None)
      }
    }
    //accumulated all lines for next read, decrement remaining reads and add to chunk
    else if (acc.size == 4) {
      //create fasta entry
      val read = fastq2FastaEntry(acc)
      //only add read if it meets minimum size threshold
      if (read.length < minsize) loadFastqChunk(iterator, List(), remaining_reads, minsize, chunk)
      else loadFastqChunk(iterator, List(), remaining_reads - 1, minsize, chunk.:+(read))
    }
    //still processing a read
    else {
      //get next line
      val current_line = iterator.next()
      //sanity check
      if (acc.isEmpty) assert(current_line.startsWith("@"), "Unexpected start of a new read: " + current_line)
      //add line to accumulator and move on
      loadFastqChunk(iterator, acc.:+(current_line), remaining_reads, minsize, chunk)
    }
  }


  /**
    * Tail-recursive method to load some maximum number of reads into a chunk from a FASTA-formatted file.
    *
    * @param iterator        FASTA iterator
    * @param chunk           Accumulating chunk of reads
    * @param remaining_reads Remaining number of reads to load
    * @param minsize         Minimum read size threshold to add
    * @param acc_name        ID of read being processed
    * @param acc_sequence    Sequence of read being processed
    * @return (List[FastaEntry], Option[String])
    */
  @tailrec final def loadFastaChunk(iterator: Iterator[String],
                                    chunk: List[FastaEntry],
                                    remaining_reads: Int,
                                    minsize: Int,
                                    acc_name: Option[String],
                                    acc_sequence: StringBuilder): (ParSeq[FastaEntry], Option[String]) = {
    //empty iterator
    if (!iterator.hasNext) {
      //no additional reads means reached end of file, check that there is a remaining read with sequence
      assert(acc_name != None, "Expected FASTA-entry at end of file")
      assert(!acc_sequence.isEmpty, "Expected corresponding sequence at end of file")
      //add last read to iterator
      val seq = acc_sequence.mkString
      (chunk.:+(new FastaEntry(acc_name.get, encodeSequence(seq), seq.size)).par, None)
    }
    //loaded max number of reads
    else if (remaining_reads == 0) {
      //check that there is a read entry in the acc
      assert(!acc_name.isEmpty, "Expected read entry after loading max number of reads: " + iterator.next())
      //check that the corresponding sequence is empty
      assert(acc_sequence.isEmpty, "Expected no corresponding sequence for read entry after loading max number of " +
        "reads: " + acc_sequence.mkString(""))
      (chunk.par, acc_name)
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
        loadFastaChunk(iterator, chunk, remaining_reads, minsize, Option(getReadName(line)), acc_sequence)
      }
      else {
        //for when there is an existing FASTA-entry and the next line is a new FASTA-entry
        if (line.startsWith(">")) {
          //sanity check
          assert(acc_sequence != None, "Expected FASTA sequence for entry " + acc_name.get + ", instead" +
            " found start of new entry: " + line)
          //only add if fasta entry meets minimum read length
          if (acc_sequence.length < minsize)
            loadFastaChunk(iterator, chunk, remaining_reads, minsize, Option(getReadName(line)), new StringBuilder)
          else {
            //get sequence
            val seq = acc_sequence.mkString
            //create a new FastaEntry object
            val fasta_entry = new FastaEntry(acc_name.get, encodeSequence(seq), seq.size)
            //add fasta entry and start accumulating new read
            loadFastaChunk(iterator, chunk.:+(fasta_entry), remaining_reads - 1, minsize, Option(getReadName(line)),
              new StringBuilder)
          }
        }
        //for when there is an existing FASTA entry and the next line is a sequence line
        else loadFastaChunk(iterator, chunk, remaining_reads, minsize, acc_name, acc_sequence.append(line))
      }
    }
  }
}

