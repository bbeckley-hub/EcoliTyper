/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#define _XOPEN_SOURCE 600
#include <limits.h>
#include <stdlib.h>
#include "hashmapkma.h"
#include "hashmapkmers.h"
#include "hashtable.h"
#include "pherror.h"

int intpos_bin(const unsigned *str1, const int str2) {
	
	int pos, upLim, downLim;
	
	upLim = *str1;
	if(upLim == 0) {
		return -1;
	}
	
	downLim = 1;
	pos = (upLim + downLim) >> 1;
	while(0 < (upLim - downLim)) {
		if(str1[pos] == str2) {
			return pos;
		} else if(str1[pos] < str2) {
			downLim = pos + 1;
		} else {
			upLim = pos - 1;
		}
		pos = (upLim + downLim) >> 1;
	}
	if(str1[pos] == str2) {
		return pos;
	}
	return -1;
}

HashTable * collect_Kmers(const HashMapKMA *templates, unsigned *Scores, long unsigned *Scores_tot, HashMap_kmers *foundKmers, Hit *hits) {
	
	int template, SU;
	long unsigned i;
	unsigned j, *value;
	short unsigned *values_s;
	HashTable_kmers *node, *node_next;
	HashTable *kmerNode, *kmerList;
	
	if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	for(i = 0; i <= foundKmers->size; ++i) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			value = hashMap_get(templates, node->key);
			if(value) {
				++hits->n;
				hits->tot += node->value;
				
				kmerNode = malloc(sizeof(HashTable));
				if(!kmerNode) {
					ERROR();
				}
				if(SU) {
					values_s = (short unsigned *) value;
					kmerNode->value = smalloc(((*values_s) + 1) * sizeof(unsigned));
					*(kmerNode->value) = *values_s;
					j = *values_s + 1;
					while(--j) {
						template = values_s[j];
						kmerNode->value[j] = template;
						Scores[template]++;
						Scores_tot[template] += node->value;
					}
				} else {
					kmerNode->value = smalloc(((*value) + 1) * sizeof(unsigned));
					*(kmerNode->value) = *value;
					j = *value + 1;
					while(--j) {
						template = value[j];
						kmerNode->value[j] = template;
						Scores[template]++;
						Scores_tot[template] += node->value;
					}
				}
				kmerNode->key = node->value;
				
				kmerNode->next = kmerList;
				kmerList = kmerNode;
			}
			free(node);
		}
	}
	free(foundKmers->table);
	
	return kmerList;
}

HashTable ** collect_Kmers_deCon(const HashMapKMA *templates, unsigned *Scores, long unsigned *Scores_tot, HashMap_kmers *foundKmers, Hit *hits, int contamination) {
	
	int template, SU;
	long unsigned i;
	unsigned j, n, *value;
	short unsigned *value_s;
	HashTable_kmers *node, *node_next;
	HashTable *kmerNode, *kmerList;
	HashTable *decon_node, *deconList;
	HashTable **Returner;
	
	if(templates->DB_size < USHRT_MAX) {
		SU = 1;
	} else {
		SU = 0;
	}
	value_s = 0;
	
	hits->n = 0;
	hits->tot = 0;
	
	kmerList = 0;
	kmerNode = 0;
	
	deconList = 0;
	decon_node = 0;
	
	Returner = smalloc(sizeof(HashTable *) << 1);
	
	for(i = 0; i <= foundKmers->size; ++i) {
		for(node = foundKmers->table[i]; node != 0; node = node_next) {
			node_next = node->next;
			if((value = hashMap_get(templates, node->key))) {
				/* check for contamination */
				++hits->n;
				hits->tot += node->value;
				
				if(SU) {
					value_s = (short unsigned *) value;
					n = *value_s;
					j = value_s[*value_s];
				} else {
					n = *value;
					j = value[*value];
				}
				
				if(j == contamination) {
					decon_node = smalloc(sizeof(HashTable));
					decon_node->value = smalloc((n + 1) * sizeof(unsigned));
					decon_node->value[0] = n;
					j = n + 1;
					if(SU) {
						while(--j) {
							decon_node->value[j] = value_s[j];
						}
					} else {
						while(--j) {
							decon_node->value[j] = value[j];
						}
					}
					decon_node->key = node->value;
					
					decon_node->next = deconList;
					deconList = decon_node;
				} else {
					kmerNode = smalloc(sizeof(HashTable));
					kmerNode->value = smalloc((n + 1) * sizeof(unsigned));
					kmerNode->value[0] = n;
					j = n + 1;
					if(SU) {
						while(--j) {
							template = value_s[j];
							kmerNode->value[j] = template;
							Scores[template]++;
							Scores_tot[template] += node->value;
						}
					} else {
						while(--j) {
							template = value[j];
							kmerNode->value[j] = template;
							Scores[template]++;
							Scores_tot[template] += node->value;
						}
					}
					
					kmerNode->key = node->value;
					
					kmerNode->next = kmerList;
					kmerList = kmerNode;
				}
			}
			free(node);
		}
	}
	free(foundKmers->table);
	
	Returner[0] = kmerList;
	Returner[1] = deconList;
	
	return Returner;
}

HashTable * withDraw_Kmers(unsigned *Scores, long unsigned *Scores_tot, HashTable *kmerList, int template, Hit *hits) {
	
	unsigned i;
	HashTable *node, *prev;
	prev = 0;
	
	if(kmerList == 0) {
		return 0;
	}
	
	node = kmerList;
	while(node != 0) {
		if(intpos_bin(node->value, template) != -1) {
			--hits->n;
			hits->tot -= node->key;
			for(i = 1; i <= node->value[0]; ++i) {
				Scores[node->value[i]]--;
				Scores_tot[node->value[i]] -= node->key;
			}
			if(prev == 0) {
				kmerList = node->next;
				free(node->value);
				free(node);
				node = kmerList;
			} else {
				prev->next = node->next;
				free(node->value);
				free(node);
				node = prev->next;
			}
			/* early stopping */
			if(Scores[template] == 0 && Scores_tot[template] == 0) {
				return kmerList;
			}
			
		} else {
			prev = node;
			node = node->next;
		}
	}
	
	if(prev != 0) {
		prev->next = 0;
	}
	
	return kmerList;
}

Hit withDraw_Contamination(unsigned *Scores, long unsigned *Scores_tot, HashTable *kmerList, HashTable *deConTable, int template, Hit hits) {
	
	unsigned i, belong;
	HashTable *node, *prev;
	HashTable *cont_node, *cont_prev;
	if(kmerList != 0) {
		prev = 0;
		node = kmerList;
		cont_prev = 0;
		cont_node = deConTable;
		while(node != 0) {
			/* check if k-mer belongs to template */
			belong = 0;
			for(i = node->value[0]; i != 0 && !belong; --i) {
				if(node->value[i] == template) {
					belong = 1;
				}
			}
			if(belong) { //withdraw score
				--hits.n;
				hits.tot -= node->key;
				for(i = node->value[0]; i != 0; --i) {
					Scores[node->value[i]]--;
					Scores_tot[node->value[i]] -= node->key;
				}
				cont_node->value = node->value;
				cont_node->key = node->key;
				cont_node->next = smalloc(sizeof(HashTable));
				cont_prev = cont_node;
				cont_node = cont_node->next;
				
				/* delete node */
				if(prev == 0) {
					kmerList = node->next;
					free(node);
					node = kmerList;
				} else {
					prev->next = node->next;
					free(node);
					node = prev->next;
				}
			} else {
				prev = node;
				node = node->next;
			}
		}
		if(cont_prev != 0) {
			free(cont_prev->next);
			cont_prev->next = 0;
		} else {
			free(deConTable);
			deConTable = 0;
		}
	} else {
		free(deConTable);
		deConTable = 0;
	}
	return hits;
}
