	/*Resort nodes and saddles in increasing order of altitude, it's sometimes used!*/
	/*This algorithm asumed both unifying lakes whre sorted as well*/
	for (i=0; i<Lake[ild].n; i++) {
	    for (j=0; j<Lake[il].n; j++) {
		if (topo[Lake[ild].row[i]][Lake[ild].col[i]] < topo[Lake[il].row[j]][Lake[il].col[j]]) {
		    /*Shift upwards the lake nodes above this one to make romm for it*/
		    for (k=Lake[il].n; k>j; k--) {
		    PRINT_DEBUG(">%d(%d,%d) %d(%d,%d)   %d   %.2f %.2f", i, Lake[ild].n, ild, j, Lake[il].n, il, k, topo[Lake[ild].row[i]][Lake[ild].col[i]], topo[Lake[il].row[j]][Lake[il].col[j]]);
		    	Lake[il].row[k] = Lake[il].row[k-1];
			Lake[il].col[k] = Lake[il].col[k-1];
		    }
		    /*Now transfer the deleted-lake node to that place*/
		    Lake[il].row[j] = Lake[ild].row[i];
		    Lake[il].col[j] = Lake[ild].col[i];
		    break;
		}
	    } 
	}
	for (i=0; i<Lake[ild].n_sd; i++) {
	    for (j=0; j<Lake[il].n_sd; j++) {
	        PRINT_DEBUG("@%d %d   %d %d   %d %d  %d %d  %d %d", i, j, ild, il, Lake[ild].n_sd, Lake[il].n_sd, Lake[ild].row_sd[i], Lake[ild].col_sd[i], Lake[14].row_sd[1], Lake[14].col_sd[1]);
	        PRINT_DEBUG("!%d %d   %d %d   %d %d  %d %d  %d %d", i, j, ild, il, Lake[ild].n_sd, Lake[il].n_sd, Lake[ild].row_sd[i], Lake[ild].col_sd[i], Lake[il].row_sd[j], Lake[il].col_sd[j]);
	        PRINT_DEBUG("#%d %d %.2f %.2f", i, j, topo[Lake[ild].row_sd[i]][Lake[ild].col_sd[i]], topo[Lake[il].row_sd[j]][Lake[il].col_sd[j]]);
		if (topo[Lake[ild].row_sd[i]][Lake[ild].col_sd[i]] < topo[Lake[il].row_sd[j]][Lake[il].col_sd[j]]) {
		    /*Shift upwards the lake nodes above this one to make romm for it*/
		    for (k=Lake[il].n_sd; k>j; k--) {
		    PRINT_DEBUG("#%d(%d,%d) %d(%d,%d)   %d   %.2f %.2f", i, Lake[ild].n, ild, j, Lake[il].n, il, k, topo[Lake[ild].row[i]][Lake[ild].col[i]], topo[Lake[il].row[j]][Lake[il].col[j]]);
		    	Lake[il].row_sd[k] = Lake[il].row_sd[k-1];
			Lake[il].col_sd[k] = Lake[il].col_sd[k-1];
		    }
		    /*Now transfer the deleted-lake node to that place*/
		    Lake[il].row_sd[j] = Lake[ild].row_sd[i];
		    Lake[il].col_sd[j] = Lake[ild].col_sd[i];
		    break;
		}
		    PRINT_DEBUG("#%d(%d,%d) %d(%d,%d)   %d   %.2f %.2f", i, Lake[ild].n, ild, j, Lake[il].n, il, k, topo[Lake[ild].row[i]][Lake[ild].col[i]], topo[Lake[il].row[j]][Lake[il].col[j]]);
	    } 
	}
