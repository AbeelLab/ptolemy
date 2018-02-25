/*
 * This work is licensed under the Creative Commons
 * Attribution-NonCommercial-NoDerivs 3.0 Unported License.
 * To view a copy of this license, visit
 * http://creativecommons.org/licenses/by-nc-nd/3.0/
 * or send a letter to Creative Commons, 444 Castro Street,
 * Suite 900, Mountain View, California, 94041, USA.
 *
 * A copy of the license is included in LICENSE.txt
 * See the License for the specific language governing permissions
 * and limitations under the License.
 *
 * Copyright 2005-2016 Thomas Abeel
 */
package atk

import java.io.{File, PrintWriter}
import java.lang.management.ManagementFactory
import java.text.{NumberFormat, SimpleDateFormat}
import java.time.LocalDateTime
import java.util.{Date, Locale, TimeZone}

object Tool extends Tool {

}

/**
  * Utility methods to create tools
  *
  * @author Thomas Abeel
  */
trait Tool extends LoggingTrait{

  val version = "There is no version information available for this tool"

  val description = "There is no description available for this tool"

  val nfP = NumberFormat.getPercentInstance(Locale.US)
  nfP.setMaximumFractionDigits(2)

  val nf = NumberFormat.getInstance(Locale.US)
  nf.setMaximumFractionDigits(2)

  val nf0 = NumberFormat.getInstance(Locale.US)
  nf0.setMaximumFractionDigits(0)

  //lazy val naturalOrdering = Ordering.comparatorToOrdering(NaturalOrderComparator.NUMERICAL_ORDER_IGNORE_CASE)
  private var logger: PrintWriter = null;

  private val timestampFormat: SimpleDateFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH-mm-ss-SSS")

  timestampFormat.setTimeZone(TimeZone.getTimeZone("UTC"))

  def timestamp(): String = {
    timestampFormat.format(new Date(System.currentTimeMillis()))

  }

  private val startTime = System.currentTimeMillis();
  private var progressCounter = 1

  def progress(reportFreq: Int) = {
    if (progressCounter % reportFreq == 0) {
      val interval = System.currentTimeMillis() - startTime
      log("--Processing: " + nf0.format(progressCounter) + "\t" + new TimeInterval(interval) + "\t" + nf.format
      (progressCounter * 1000L / (interval + .1)) + " units/s\t" + LocalDateTime.now())
    }
    progressCounter += 1
  }

  def init(location: String = null, config: AnyRef = null): Unit = {
    if (location == null){
      logger = new PrintWriter(System.out)
      //logger.print(generatorInfo(config) + "\n")
    }
    else
      logger = new NixWriter(location,config)

  }

  def log(str: Any) = {
    if (logger == null) {
      init()

    }
    logger.println(str)
    logger.flush()
  }

  def finish(out: PrintWriter = logger) = {
    print("## This analysis finished " + new Date() + "\n")
    if (out != null) {
      out.print("## This analysis finished " + new Date() + "\n")
      out.print("## Run time: " + new TimeInterval(System.currentTimeMillis() - startTime) + "\n")
      out.close()
    }
  }
  private def classFileName() = { Thread.currentThread().getStackTrace().takeRight(3).map(_.getFileName()).mkString(";") };

  private def classInfo() = { Thread.currentThread().getStackTrace().takeRight(3).map(_.getClassName()).mkString(";") };

  private def executeEnvironment() = { this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath() }

  private def getDeclaredFields(cc: AnyRef) = {
    val m = (Map[String, Any]() /: cc.getClass.getDeclaredFields) { (a, f) =>
      f.setAccessible(true)
      a + (f.getName -> f.get(cc))
    }
    m.toList.sortBy(_._1)
  }

  /* Added for legacy purposes */
  def generatorInfo(): String = { generatorInfo(null) };

  def generatorInfo(config: AnyRef = null): String = {
    "\n# Generated with    " + classInfo() + "\n" +
      "# Source code in    " + classFileName() + "\n" +
      "# Binary in         " + executeEnvironment + "\n" +
      "# Date              " + new Date() + "\n" +
      "# Working directory " + new File(".").getAbsolutePath() + "\n#\n" +
      "# Machine configuration summary: \n#\t " + "Current date and time: " + new Date() + "\n#\t " + "Number of processors: " + Integer.toString(Runtime.getRuntime().availableProcessors()) + "\n#\t " + "Free memory :" + Runtime.getRuntime().freeMemory() + "\n#\t " + "Max memory: " + Runtime.getRuntime().maxMemory() + "\n#\t " + "Total JVM: " + Runtime.getRuntime().totalMemory() + "\n#\t " + "OS: " + ManagementFactory.getOperatingSystemMXBean().getName() + " " + ManagementFactory.getOperatingSystemMXBean().getVersion() + "\n#\t " + "Architecture: " + ManagementFactory.getOperatingSystemMXBean().getArch() + "\n#\t " + "JVM version: " + System.getProperty("java.version") + "\n"

  }
  /* Elapsed time function */
  def time[R](block: => R): R = {
    val t0 = System.currentTimeMillis()
    val result = block // call-by-name
    val t1 = System.currentTimeMillis()
    println("Elapsed time: " + (t1 - t0) + "ms")
    result
  }

}