LIB:=lib/tailcompression.cpp lib/nloglike.cpp lib/spr.cpp lib/multinomial.cpp lib/robust_multi.cpp lib/dirichlet.cpp lib/equilibrium.cpp
INC:=include/*.h include/*.hpp

ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: build/debug/wt build/release/wt

build/debug/wt: build/debug wt.tpl 
	@admb -g wt.tpl 
	@$(COPY) wt build/debug/
	@$(RM) wt

build/release/wt: build/release wt.tpl 
	@admb -g wt.tpl 
	@$(COPY) wt build/release/
	@$(RM) wt

build/debug:
	@echo creating debug directory
	@mkdir -p build/debug

build/release:
	@echo creating release directory
	@mkdir -p build/release

clean:
	@$(RM) build
	@$(RM) wt.cpp
	@$(RM) wt.htp
	@$(RM) wt.obj
