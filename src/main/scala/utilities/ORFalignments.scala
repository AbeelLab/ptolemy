package utilities

/**
  * Author: Alex N. Salazar
  * Created on 14-2-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait ORFalignments {

  /**
    * Case class for pairwise ORF alignment
    * @param query
    * @param ref
    * @param sigma
    * @param gamma
    */
  case class ORFalignment(query: Int, ref: Int, sigma: Double, gamma: Double)

  /**
    * Function to determine
    * @return
    */
  def computeSigma: Array[String] => Double = alignment =>alignment(9).toDouble/alignment(1).toDouble


  /**
    * Function to determine Gamma (min seq length / max seq length)
    * @return
    */
  def computeGamma: Array[String] => Double = alignment => {
    //alignment(10).toDouble/alignment(1).toDouble
    val (a,b) = (alignment(1).toDouble, alignment(6).toDouble)
    List(a,b).min / List(a,b).max
  }

  /**
    * Function to convert pairwise alignment to ORFalignment object
    * @return
    */
  def toORFalignment: Array[String] => ORFalignment = alignment =>
    new ORFalignment(alignment(0).toInt, alignment(5).toInt, computeSigma(alignment), computeGamma(alignment))

}
