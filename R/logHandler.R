#' class for handling logs/warnings/errorss
#'
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @export
#' @import stringr XML
logHandler<-function() {
  sendLog<-function( msg , type="MSG") {
    if( type == "MSG") cat("\n",msg)
    else stop( msg )
  }
  costructor<-function() {
  }
  costructor();
  return(
    list(
      "sendLog"=sendLog
    )
  )
}
