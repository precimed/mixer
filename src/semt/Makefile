.ONESHELL:
SHELL=/bin/bash

EMPTY :=
SPACE := ${EMPTY} ${EMPTY}

BUILD_DIR := ./Debug
GCOV_DIR := ./gcov
DOXY_DIR := ./doc

DIRS := semt/ examples/ loki/
SEMT_FILES = $(foreach dir,$(DIRS),$(wildcard $(dir)/*))
VPATH := ${DIRS}

#
# compiler settings
#
CC := g++
CFLAGS := -Wall -std=c++0x -g3 -fmessage-length=0 -ftemplate-depth-200
#-DSEMT_DISABLE_PRINT #-DSEMT_DISABLE_ALL_SIMPLIFICATIONS
CBFLAGS := -O0
INCLUDES := -I ./
LIBS :=

#
# predicates
#
ADDITIONAL_COMMANDS = @set -e;

ifneq (,$(findstring release,$(MAKECMDGOALS)))
	BUILD_DIR := ./Release
	CBFLAGS := -O2 -DNDEBUG
	ADDITIONAL_COMMANDS += strip --strip-unneeded $@
endif

ifneq (,$(findstring pedantic,$(MAKECMDGOALS)))
	CBFLAGS += -pedantic
endif

ifneq (,$(findstring gprof,$(MAKECMDGOALS)))
	CBFLAGS += -pg
	ADDITIONAL_COMMANDS += ${BUILD_DIR}/${TARGET} $(ARGS)|| true;\
		mkdir -p ./gprof;\
		gprof ${BUILD_DIR}/${TARGET} > ./gprof/$(subst ${SPACE},_,${TARGET}_$(ARGS));\
		rm ${BUILD_DIR}/*.o;
endif

ifneq (,$(findstring gcov,$(MAKECMDGOALS)))
	CBFLAGS += -fprofile-arcs -ftest-coverage
	LIBS += -lgcov
	ADDITIONAL_COMMANDS += ${BUILD_DIR}/${TARGET} $(ARGS);\
		gcov -o ${BUILD_DIR} ${SRC_FILE} > /dev/null;\
		mkdir -p ${GCOV_DIR};\
		mv *.gcov ${GCOV_DIR};\
		rm ${BUILD_DIR}/*.o;
endif

#
# objects
#
SEMTOBJS := ${BUILD_DIR}/VectorExpr.o
semt_check_OBJS := ${SEMTOBJS}
semt_speed_OBJS := ${BUILD_DIR}/semt_speed_func.o ${SEMTOBJS}
semt_newton_OBJS := ${BUILD_DIR}/semt_func_impl.o ${SEMTOBJS}
semt_examples_OBJS := ${SEMTOBJS}

BINARIES := semt_func_impl \
		semt_check \
		semt_newton \
		semt_speed_func \
		semt_speed \
		semt_examples

TARGET := $(filter ${BINARIES},$(strip $(MAKECMDGOALS)))
SRC_FILE := $(foreach dir,$(DIRS),$(wildcard $(dir)/${TARGET}.cpp))

#
# dependencies
#
DEPS := $(patsubst %,${BUILD_DIR}/%.d,${BINARIES})
-include $(DEPS)

${BUILD_DIR}/%.d: %.cpp | $(BUILD_DIR)
	@set -e; rm -f $@; \
	$(CC) ${INCLUDES} ${CFLAGS} ${CBFLAGS} -MM -MF"$@.tmp" $<; \
	sed 's,\($*\)\.o[ :]*,$(BUILD_DIR)/\1.o $@ : ,g' < $@.tmp > $@; \
	rm -f $@.tmp;

${BUILD_DIR}/%.o: %.cpp | $(BUILD_DIR)
	${CC} ${INCLUDES} -c ${CFLAGS} ${CBFLAGS} -o"$@" $<

$(BUILD_DIR):
	@mkdir -p $@

#
# targets
#
.PHONY: all build debug release gcov gprof touch_target \
		clean distclean \
		doc ${BINARIES}
.PRECIOUS: ${BUILD_DIR}/%.o

all build: semt_check semt_examples semt_speed semt_newton

debug release:
	@echo "Using $@ configuration"

gcov gprof: touch_target ${TARGET}
	@echo "done $@"

doc: ${SEMT_FILES} doxy_pages Doxyfile
	doxygen

clean:
	@rm -rf ./Debug ./Release
	@rm -rf ${GCOV_DIR}
	@find . -name "*.gc*" -delete
	@find . -name "*.o" -delete
	@find . -name "*.out" -delete
	@echo "cleaned project folder"

distclean: clean
	@rm -rf ${DOXY_DIR}
	@find . -name "*~" -delete

${BINARIES} : % : ${BUILD_DIR}/%

.SECONDEXPANSION:
$(BUILD_DIR)/%: $(BUILD_DIR)/%.o $$($$(*F)_OBJS)
	@echo "linking $^"
	@${CC} ${INCLUDES} ${CFLAGS} ${CBFLAGS} -o $@ $^
	$(ADDITIONAL_COMMANDS)

touch_target:
	@touch ${SRC_FILE}
	@for file in $(${TARGET}_OBJS); do \
		[ ! -f $$file ] || rm $$file; \
	done

%:
	@true;
