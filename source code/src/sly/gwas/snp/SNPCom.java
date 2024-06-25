package sly.gwas.snp;

import sly.gwas.mem.Mem;

import java.util.Arrays;

public class SNPCom implements Comparable<SNPCom>{
    int[] x=null;
    public SNPCom(int[] x){
        this.x=x;
    }

    @Override
    public int compareTo(SNPCom o) {
        if (this.x.length < o.x.length) {
            return -1;
        } else if (this.x.length > o.x.length) {
            return 1;
        } else {
            for (int i = 0; i < x.length; i++) {
                if (x[i] < o.x[i]) {
                    return -1;
                } else if (x[i] > o.x[i]) {
                    return 1;
                }
            }
            return 0;
        }
    }

    @Override
    public String toString() {
        return Arrays.toString(x);
    }

    public boolean contains(int[] e){
        int same=0;
        for(int i=0,j=0;i<x.length&&j<e.length;){
            if(x[i]==e[j]){
                same++;
                i++;
                j++;
            }
            else if(x[i]<e[j]){
                i++;
            }
            else{
                j++;
            }
        }
        return same==e.length;
    }
}
