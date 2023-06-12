
gpt4.gcc: gpt4.c
	gcc -I ~/src -I ~/src/arb -c gpt4.c 2>&1 | tee gpt4.gcc
