##
## post-10-mclapply.hack.R
##
## Nathan VanHoudnos
## nathanvan AT northwestern FULL STOP edu
## July 14, 2014
## Last Edit:  August 26, 2014
##
## A script to implement a hackish version of
## parallel:mclapply() on Windows machines.
## On Linux or Mac, the script has no effect
## beyond loading the parallel library.

## Define the hack
mclapply.hack <- function(..., mc.cores=NULL) {
    ## Create a cluster
    if( is.null(mc.cores) ) {
        size.of.list <- length(list(...)[[1]])
        mc.cores <- min(size.of.list, detectCores())
    }
    ## N.B. setting outfile to blank redirects output to
    ##      the master console, as is the default with
    ##      mclapply() on Linux / Mac
    cl <- makeCluster( mc.cores, outfile="" )

    ## Find out the names of the loaded packages
    loaded.package.names <- c(
        ## Base packages
        sessionInfo()$basePkgs,
        ## Additional packages
        names( sessionInfo()$otherPkgs ))

    tryCatch( {

        ## Copy over all of the objects within scope to
        ## all clusters.
        this.env <- environment()
        while( identical( this.env, globalenv() ) == FALSE ) {
            clusterExport(cl,
                          ls(all.names=TRUE, env=this.env),
                          envir=this.env)
            this.env <- parent.env(environment())
        }
        clusterExport(cl,
                      ls(all.names=TRUE, env=globalenv()),
                      envir=globalenv())

        ## Load the libraries on all the clusters
        ## N.B. length(cl) returns the number of clusters
        parLapply( cl, 1:length(cl), function(xx){
            lapply(loaded.package.names, function(yy) {
                require(yy , character.only=TRUE)})
        })

        ## Run the lapply in parallel
        return( parLapply( cl, ...) )
    }, finally = {
        ## Stop the cluster
        stopCluster(cl)
    })
}

## Warn the user if they are using Windows
if( Sys.info()[['sysname']] == 'Windows' ){
    message(paste(
        "\n",
        "   *** Microsoft Windows detected ***\n",
        "   \n",
        "   For technical reasons, the MS Windows version of mclapply()\n",
        "   is implemented as a serial function instead of a parallel\n",
        "   function.",
        "   \n\n",
        "   As a quick hack, we replace this serial version of mclapply()\n",
        "   with a wrapper to parLapply() for this R session. Please see\n\n",
        "     http://www.stat.cmu.edu/~nmv/2014/07/14/implementing-mclapply-on-windows \n\n",
        "   for details.\n\n"))
}

## If the OS is Windows, set mclapply to the
## the hackish version. Otherwise, leave the
## definition alone.
mclapply <- switch( Sys.info()[['sysname']],
                    Windows = {mclapply.hack},
                    Linux   = {mclapply},
                    Darwin  = {mclapply})

## end post-10-mclapply.hack.R
