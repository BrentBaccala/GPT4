
# Programming with GPT4

Hi!

For years I've been thinking that in the far future, many years after my demise,
computers will have learned to understand English, and we'll be able to program
them like they do in Star Trek.

I assumed that I'd never live to see this, and people would look back on my time
as the dark ages of computing, when we had to learn arcane languages like C++
and Java.  Only a small percentage of the population could really program a
computer, and those people slaved for years to produce bug-riddled code
that did enough to be useful.

Then along came Large Language Models!

I'd been hearing about the transformer architecture for years, but I'd never
used one.  I knew the A.I. guys really liked it, but A.I. isn't my specialty.
I work in differential algebra and computational algebraic geometry,
stuff like that.

Now that I've started using GPT-4, I see that this must have been low
hanging fruit for at least the last five years.  The A.I. experts knew
that they had a very nice neural network architecture for processing
human language text.  An obvious thing to try was to build a really
big one and train it on the entire Internet.  We're ten years into
the cloud revolution, so if you've got enough money, Microsoft,
Amazon, Google, and a dozen other major players have got data centers
all over the world packed with high end servers and GPUs that you can rent over
the Internet.

The results are nothing short of amazing.

and remember, I write software!

So here we go...

### awk script

Here's a simple GPT-4 program that produces an awk script to add up my timesheet:

```
My work log contains some lines that start will billable hours.
All such lines start with optional whitespace, then a decimal number of hours, which may start with a decimal, then the word 'hr'.
Write an awk script that adds up all billable hours in my work log.
Make sure the script starts with a shebang command.
Output the script without commentary.
```

Here's GPT-4's output:

```
#!/usr/bin/awk -f

BEGIN { total = 0 }
/^[ \t]*([0-9]+(\.[0-9]+)?|[.][0-9]+) +hr/ { total += $1 }
END { printf "Total billable hours: %.2f\n", total }
```

### web page scrapper

I've been using GPT-4 on poe.com, and wrote this GPT-4 program to automate the process.

It was also a nice programming exercise.

- [GPT-4 prompt](automate)
- [GPT-4 output](automate.out)

It doesn't completely work, because I can't figure how to get the Copy button to scroll
into view so that it can be clicked on, but I can't figure it out myself, so how could the machine?

### Grobner basis

Here's my most complicated GPT-4 program to date.

It computes a non-reduced Grobner basis using the naive Buchberger's algorithm
as documented on [this wikipedia page](wikipedia.org/Buchberger's_Algorithm).

There is a single test case, taken from [Bernd Sturmfel's article](https://math.berkeley.edu/~bernd/what-is.pdf)

The input prompts are pre-processed using `m4`, so GPT-4 never sees an
`include` line, but rather gets the common file included.

- [GPT-4 "include" file](common)
- [GPT-4 input5](input5)
- [GPT-4 C code output](input5.c)
- [GPT-4 input1](input1)
- [GPT-4 C code output](input1.c)
- [GPT-4 input2](input2)
- [GPT-4 C code output](input2.c)
- [GPT-4 input3](input3)
- [GPT-4 C code output](input3.c)
- [GPT-4 input4](input4)
- [GPT-4 C code output](input4.c)
- [Makefile](Makefile)
- [program output](gpt4.out)

## Conclusion

I'm learning to program in GPT-4!  It's close enough to English that I have
hopes for it to be fairly easily portable to GPT-5.  I also think that
it's the programming style of the future.
