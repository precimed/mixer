#!/bin/bash
OUTDIR=./Debug
SRCDIR=./examples
EXE=semt_main
OBJS="semt_speed_func semt_speed"
TEMP_FILE=$OUTDIR/evil_temp_file

OPTS="0 1"

CFLAGS="-fmessage-length=0 -MMD -MP -Wall -std=c++0x -g3 -ftemplate-depth-200 -DNDEBUG"

mkdir -p $OUTDIR

for OPT in $OPTS;
do 
    echo "Optimization level $OPT:"
    TOTAL_TIME=0
    COMPILE_TIME=0
    RUN_TIME=0
    
    echo -ne "\tcompiling vector... "
    /usr/bin/time -f "%U" -o $TEMP_FILE \
    	g++  -I ./ -o"$OUTDIR/semt_vector$OPT.o" -c $CFLAGS -O$OPT -MF"$OUTDIR/semt_vector$OPT.d" -MT"$OUTDIR/semt_vector$OPT.d" semt/VectorExpr.cpp &> /dev/null;
    COMPILE_TIME=$(cat $TEMP_FILE)
    echo "done. ($COMPILE_TIME)"
    
    for OBJ in $OBJS;
    do
    	echo -ne "\tcompiling $OBJ... "
    	/usr/bin/time -f "%U" -o $TEMP_FILE \
    		g++  -I ./ -o"$OUTDIR/$OBJ$OPT.o" -c $CFLAGS -O$OPT -MF"$OUTDIR/$OBJ$OPT.d" -MT"$OUTDIR/$OBJ$OPT.d" $SRCDIR/$OBJ.cpp &> /dev/null;
    	RUN_TIME=$(cat $TEMP_FILE)
    	COMPILE_TIME=$(echo "$COMPILE_TIME + $RUN_TIME" | bc)
    	echo "done. ($RUN_TIME)"
    done
    
    echo -ne "\tlinking... "
    /usr/bin/time -f "%U" -o $TEMP_FILE \
    	g++ -o $OUTDIR/$EXE$OPT $CFLAGS -O$OPT $OUTDIR/semt_vector$OPT.o $OUTDIR/semt_speed_func$OPT.o $OUTDIR/semt_speed$OPT.o &> /dev/null;
    RUN_TIME=$(cat $TEMP_FILE)
    COMPILE_TIME=$(echo "$COMPILE_TIME + $RUN_TIME" | bc)
    echo "done. ($RUN_TIME)"
    
    echo -ne "\trunning... "
    /usr/bin/time -f "%U" -o $TEMP_FILE $OUTDIR/$EXE$OPT > /dev/null;
    RUN_TIME=$(cat $TEMP_FILE)
    echo "done. ($RUN_TIME)"
    
    TOTAL_TIME=$(echo "$COMPILE_TIME + $RUN_TIME" | bc)
    echo -e "overall: $TOTAL_TIME s (compile: $COMPILE_TIME s, run: $RUN_TIME s)\r\n\r\n"
done
rm $TEMP_FILE
