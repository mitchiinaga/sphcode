CXX = g++

# C++11
CXXFLAGS = -std=gnu++14

# Optimization
CXXFLAGS += -O3 -funroll-loops -ffast-math

# OpenMP
CXXFLAGS += -fopenmp

# use for debug
CXXFLAGS += -Wall
CXXFLAGS += -Wno-sign-compare
CXXFLAGS += -Wno-maybe-uninitialized
#CXXFLAGS += -g =D_DEBUG

SRCDIR = src
BUILDDIR = obj
INC = -Iinclude

SRCEXT = cpp
DEPEXT = d
OBJEXT = o

TARGET = sph

# DO NOT EDIT BELOW THIS LINE
#-------------------------------------------------------

sources = $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
objects = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(subst $(SRCEXT),$(OBJEXT),$(sources)))
dependencies = $(subst .$(OBJEXT),.$(DEPEXT),$(objects))

# Defauilt Make
all: $(TARGET)

# Remake
remake: clean all

# ディレクトリ生成
directories:
	@mkdir -p $(BUILDDIR)

# 中間生成物のためのディレクトリを削除
clean:
	@$(RM) -f $(BUILDDIR)/* $(TARGET)

# 自動抽出した.dファイルを読み込む
-include $(dependencies)

# オブジェクトファイルをリンクしてバイナリを生成
$(TARGET): $(objects)
	$(CXX) -o $(TARGET) $(CXXFLAGS) $^ $(FLAGS)

# ソースファイルのコンパイルしてオブジェクトファイルを生成
# また、ソースファイルの依存関係を自動抽出して.dファイルに保存
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INC) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) $(INC) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

# Non-File Targets
.PHONY: all remake clean
