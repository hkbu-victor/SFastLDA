all:
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"super_fast_final1_d.d" -MT"super_fast_final1_d.d" -o "super_fast_final1_d.o" "super_fast_final1_d.cpp"
	g++  -o "super_fast_final1_d"  super_fast_final1_d.o 
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"super_fast_final1_w.d" -MT"super_fast_final1_w.d" -o "super_fast_final1_w.o" "super_fast_final1_w.cpp"
	g++  -o "super_fast_final1_w"  super_fast_final1_w.o 
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"doc_uci_to_word_uci.d" -MT"doc_uci_to_word_uci.d" -o "doc_uci_to_word_uci.o" "doc_uci_to_word_uci.cpp"
	g++  -o "doc_uci_to_word_uci"  doc_uci_to_word_uci.o
