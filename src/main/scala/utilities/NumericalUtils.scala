package utilities

/**
  * Author: Alex N. Salazar
  * Created on 18-8-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait NumericalUtils {

  /**
    * Compute max value given two Ints
    * @return Int
    */
  def max: (Int, Int) => Int = (x, y) => if(x > y) x else y

  /**
    * Compute min value given two Ints
    * @return
    */
  def min: (Int,Int) => Int = (x, y) => if(x < y) x else y

  /**
    * Compute absolute value of a given an Int without using java's Math package
    * @return Int
    */
  def abs: Int => Int = x => if (x > 0) x else x * -1

  /**
    * Function to compute the median value given a list of ints
    *
    * @return Double (median)
    */
  def computeMedian: Iterable[Int] => Int = values => {
    def getMedian: Int => Int = midpoint => {
      val sorted = values.toSeq
      //check whether its an odd or even list
      sorted.size % 2 match {
        //if even
        case 0 => {
          //compute median value based on the average of the midpoint index and midpoint index + 1
          (sorted(midpoint) + sorted(midpoint + 1)) / 2
        }
        //if odd
        case _ => sorted((midpoint - 0.5).toInt)
      }
    }
    //get mid point index
    val mid_index = (values.size / 2.0) - 1
    //get median value
    val median = getMedian(mid_index.toInt)
    median
  }


}
