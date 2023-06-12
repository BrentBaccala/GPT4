
gpt4 gpt4.gcc: gpt4.c
	gcc -I ~/src -I ~/src/arb -g -o gpt4 gpt4.c -lflint 2>&1 | tee gpt4.gcc
