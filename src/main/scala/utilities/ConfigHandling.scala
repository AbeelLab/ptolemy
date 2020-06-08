package utilities

import java.io.{File, FileInputStream}
import java.util.Calendar
import java.util.zip.GZIPInputStream
import java.io.{File, PrintWriter}

import utilities.FileHandling.{verifyFile, timeStamp}

import scala.io.Source
import scala.util.Try

object ConfigHandling {

	case class fullConfig(
                          brh: File = null,
                          canonicalQuiver: File = null,
                          database: File = null,
                          genomesFile: File = null,
                          hlgg: File = null,
                          outputDir: File = null,
                          reads: File = null,
                          reference: File = null,
                          syntenicAnchors: File = null,
                          
                          chunkSize: Int = 30000,
                          f: Int = 10, // backcompatibility with SyntenicAnchorsC
                          flankingWindow: Int = 10,
                          kmerSize: Int = 11,
                          maxThreads: Int = 1,
                          minCoverage: Int = 5,
                          minGeneLength: Int = -1,
                          minHits: Int = 1,
                          minInterSize: Int = 11,
                          minimizerWindow: Int = 5,
                          minReadLength: Int = 500,
                          repeatRadius: Int = 10,
                          windowSize: Int = 3,
                          
                          gamma: Double = 0.75,
                          lengthProportion: Double = 0.3,
                          maxError: Double = 0.20,
                          proportion: Double = 0.6,
                          splitOverlaps: Double = 0.15,
                          syntenicFraction: Double = 0.5,
                          
                          debug: Boolean = false,
                          dump: Boolean = false,
                          isCircular: Boolean = false,
                          mergeCC: Boolean = false,
                          msa: Boolean = false,
                          showWarnings: Boolean = false,
                          traverseNodes: Boolean = false,
                          useGene: Boolean = true,
                          verbose: Boolean = false,
                          
                          prefix: String = null,
                     	 )

	// function to get the parameters of the (module) Config case class
	def mapConfigParams(cc: AnyRef) : Map[String, String] =
	cc.getClass.getDeclaredFields.foldLeft(Map.empty[String, String]) 
	{ 
		(a, f) => f.setAccessible(true) 
		a + (f.getName -> (f.get(cc)+""))
	}

	// function to save the parameters of the (module) Config case class
	def saveMapConfigParams(file:java.io.File, params: Map[String, String]) : Unit =
	{
		// join the map entries into strings stored in array
		val str = for ( (k,v) <- params) yield { (s"$k\t$v") }
		// create the file to store the configuration parameters
		val pw = new java.io.PrintWriter(file)
		// write each array entry into a line and close the file
		try pw.write(str.mkString("\n")) finally pw.close()
	}

	// function to read the parameters of the configuration file
	def parseMapConfigParams(file: String) : Map[String, String] =
	{
		// loop over the file lines ad get an array of them
		val list = for {line <- Source.fromFile(file).getLines}
			yield (line.split("\t")(0) -> line.split("\t")(1))
		list.toMap
	}

	// generate new config object with new parameters
	def generateFullConfig(params: Map[String, String]) = 
	{
		fullConfig(	
                  new File(params.get("brh").get), 
                  new File(params.get("canonicalQuiver").get),
                  new File(params.get("database").get),
                  new File(params.get("genomesFile").get), 
                  new File(params.get("hlgg").get), 
                  new File(params.get("outputDir").get),
                  new File(params.get("reads").get), 
                  new File(params.get("reference").get), 
                  new File(params.get("syntenicAnchors").get), 
                  
                  params.get("chunkSize").get.toInt,
                  params.get("f").get.toInt,
                  params.get("flankingWindow").get.toInt,
                  params.get("kmerSize").get.toInt,
                  params.get("maxThreads").get.toInt,
                  params.get("minCoverage").get.toInt,
                  params.get("minGeneLength").get.toInt,
                  params.get("minHits").get.toInt,
                  params.get("minInterSize").get.toInt,
                  params.get("minimizerWindow").get.toInt,
                  params.get("minReadLength").get.toInt,
                  params.get("repeatRadius").get.toInt,
                  params.get("windowSize").get.toInt,
                  
                  params.get("gamma").get.toDouble,
                  params.get("lengthProportion").get.toDouble,
                  params.get("maxError").get.toDouble,
                  params.get("proportion").get.toDouble,
                  params.get("splitOverlaps").get.toDouble,
                  params.get("syntenicFraction").get.toDouble,
                  
                  params.get("debug").get.toBoolean,
                  params.get("dump").get.toBoolean,
                  params.get("isCircular").get.toBoolean,
                  params.get("mergeCC").get.toBoolean,
                  params.get("msa").get.toBoolean,
                  params.get("showWarnings").get.toBoolean,
                  params.get("traverseNodes").get.toBoolean,
                  params.get("useGene").get.toBoolean,
                  params.get("verbose").get.toBoolean,
                  
                  params.get("prefix").get.toString
            )
	} 

  // two different config Map[String,String] are compared currentMap and savedMap
  // having as output another one called finalMap
	def compareMappings(currentMap: Map[String, String], savedMap: Map[String, String], module: String) : Map[String, String] =
	{
    // both maps should have the same keys, raise error if not
    val diffCS = (currentMap.keySet -- savedMap.keySet)
    val diffSC = (savedMap.keySet -- currentMap.keySet)
    val diff = diffCS ++ diffSC
    assert(diff.isEmpty, "Invalid stored config.tsv file")

    // define map of flags and the first module where they appear
    val flagModule = Map(
                          "isCircular" -> "extract",
                          "genomesFile" -> "extract",
                          "minInterSize" -> "extract",
                          "repeatRadius" -> "extract",
                          "splitOverlaps" -> "extract",
                          "showWarnings" -> "extract",
                          "useGene" -> "extract",
                          "brh" -> "syntenic-anchors", 
                          "flankingWindow" -> "syntenic-anchors", 
                          "minimizerWindow" -> "syntenic-anchors", 
                          "gamma" -> "syntenic-anchors", 
                          "syntenicFraction" -> "syntenic-anchors", 
                          "mergeCC" -> "syntenic-anchors", 
                          "kmerSize" -> "syntenic-anchors", 
                          "syntenicAnchors" -> "canonical-quiver",
                          "msa" -> "canonical-quiver",
                          "canonicalQuiver" -> "index-graph",
                          "windowSize" -> "index-graph",
                          "reads" -> "align-reads",
                          "chunkSize" -> "align-reads",
                          "minCoverage" -> "align-reads",
                          "minGeneLength" -> "align-reads",
                          "minHits" -> "align-reads",
                          "minReadLength" -> "align-reads",
                          "lengthProportion" -> "align-reads",
                          "maxError" -> "align-reads",
                          "proportion" -> "align-reads",
                          "prefix" -> "align-reads",
                          "hlgg" -> "reference-graph",
                          "reference" -> "reference-graph", 
                          "traverseNodes" -> "variant-caller", 
                        )

    // define the output map
    val finalMap = scala.collection.mutable.Map[String,String]() 

    // iterate over the keys of the map
    for ((k,v) <- currentMap) {
      // if nonmatching parameters between maps and these parameters cannot be 
      // modified between modules
      if (currentMap(k) != savedMap(k) && flagModule.contains(k)) {
        // if it is the module where the immutable parameter is first used it
        // should be chosen from input
        if (flagModule(k) == module) {
          finalMap(k) = currentMap(k)
        // otherwise savedMap values should superseed the values of currentMap
        } else {
          if (currentMap(k) != "null") {
            println(timeStamp+"Option "+k+" with value "+currentMap(k)+
                  " conflicts with previously used value "+savedMap(k)+
                  ", taking latter parameter")}
          finalMap(k) = savedMap(k)
        }
      } 
      // otherwise, just keep the currentMap (which is either the default value
      // or the user specified) or the same as savedMap if no mismatch 
      else {
        finalMap(k) = currentMap(k)
      }
    }
    return finalMap.toMap
	}

	def parameterManager(cc: AnyRef, module: String) : ConfigHandling.fullConfig =
	{
    // convert currrent module configuration into a map
	  val currentParams = mapConfigParams(cc)
    // get the path of the database and the configuration file altogether
    val configPath = new File(currentParams.get("database").getOrElse("./").toString,
                            "/config.tsv")

    // if not Extract module, config should exist, so load previously saved file
    val finalParams = if (module != "extract") 
      {
         verifyFile(configPath, "The config file cannot be found ")
	       val savedParams = parseMapConfigParams(configPath.toString)
         // compare the values of the maps and solve conflicts
         compareMappings(currentParams, savedParams, module)
      // otherwise, no saved configuration parameters
      } else { currentParams }

    // save and return the config object
	  saveMapConfigParams(configPath, finalParams)
    return generateFullConfig(finalParams)
	}

}
