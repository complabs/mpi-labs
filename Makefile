
all: 
	@echo -e "\nlab 1:"
	@$(MAKE) --no-print-directory -C lab1 all
	@echo -e "\nlab 2:"
	@$(MAKE) --no-print-directory -C lab2 all
	@echo -e "\nlab 3:"
	@$(MAKE) --no-print-directory -C lab3 all
	@echo

clean: 
	@echo -e "\nlab 1:"
	@$(MAKE) --no-print-directory -C lab1 clean
	@echo -e "\nlab 2:"
	@$(MAKE) --no-print-directory -C lab2 clean
	@echo -e "\nlab 3:"
	@$(MAKE) --no-print-directory -C lab3 clean
	@echo

run: 
	@( cd lab1; make -q; sbatch test-ex.sh )
	@( cd lab2; make -q; sbatch test-ex.sh )
	@( cd lab3; make -q; sbatch test-ex.sh )
	@qtop
	@cat lab1/test-ex.out
	@cat lab2/test-ex.out
	@cat lab3/test-ex.out
