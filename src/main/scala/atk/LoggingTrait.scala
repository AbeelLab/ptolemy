package atk

import java.util.logging.{Level, LogManager}

/**
  * Author: Alex N. Salazar
  * Created on 17-11-2017
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
trait LoggingTrait {

  setDebugLevel(Level.INFO)

  def setDebugLevel(newLvl: Level) {
    val rootLogger = LogManager.getLogManager().getLogger("");
    val handlers = rootLogger.getHandlers();
    rootLogger.setLevel(newLvl);
    for (h <- handlers) {

      h.setLevel(newLvl);
    }
  }
}
