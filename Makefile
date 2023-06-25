
INPUTS=input1.o input1a.o input1b.o input1c.o input2.o input3.o input4.o input5.o input6.o

gpt4::
	rm -f $(INPUTS)

gpt4:: $(INPUTS)
	gcc -I ~/src -I ~/src/arb -g -o gpt4 $^ -lflint

input%.o: input%.c
	gcc -I ~/src -I ~/src/arb -g -c $<
