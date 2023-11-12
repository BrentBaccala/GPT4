
# Makefile to call the GPT-4 API to convert GPT prompts into C code.

# We need to use /bin/bash to use pipefail, which is a bash option that causes a pipeline
# to fail if any command in the pipeline fails.  The command in question is "jq", which
# takes the --exit-status option to fail if its output value is null.  That's what I
# want to do if there's some kind of problem (like authentication) with the API call.

SHELL := /bin/bash

URL := https://api.openai.com/v1/chat/completions

INPUTS=input1.o input1a.o input1b.o input1c.o input2.o input3.o input4.o input5.o input6.o

gpt4:: $(INPUTS)
	gcc -I ~/src -I ~/src/arb -g -o gpt4 $^ -lflint

input%.o: input%.c
	gcc -I ~/src -I ~/src/arb -g -c $<

input%.c: input%.txt common
ifndef OPENAI_API_KEY
	$(error OPENAI_API_KEY is undefined)
endif
	set -o pipefail && m4 $< | jq -n --rawfile content /dev/stdin \
	  '{ "model": "gpt-4", "messages": [{"role": "user", "content": $$content}], "temperature": 0  }' \
	  | curl --silent $(URL) -H "Content-Type: application/json" -H "Authorization: Bearer $$OPENAI_API_KEY" -d @- \
	  | jq --exit-status -r .choices[0].message.content \
	  | sed -n '/^```c$$/,/^```$$/{//!p}' > $@
