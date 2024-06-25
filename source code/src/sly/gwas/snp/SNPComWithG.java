package sly.gwas.snp;

import sly.gwas.mem.Mem;

import java.util.Arrays;

public class SNPComWithG implements Comparable<SNPComWithG>{
    int[] x=null;
    double g=0;
    public SNPComWithG(int[] x, double g){
        this.x=x;
        this.g=g;
    }
    public final void setG(double g){
        this.g=g;
    }
    @Override
    public int compareTo(SNPComWithG o) {
        if(g<o.g){
            return 1;
        }
        else if(g>o.g){
            return -1;
        }
        else {
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
    }

    @Override
    public String toString() {
        return "SNPCom{" +
                "x=" + Arrays.toString(x) +
                ", g=" + g +
                '}';
    }

    public final String getOString(Mem mem) {
        StringBuffer r=new StringBuffer();
        for(int i=0;i<x.length;i++){
            r.append(mem.getName(x[i]));
            r.append(",");
        }
        r.append(g);
        r.append("\n");
        return r.toString();
    }

    public final double getG(){
        return this.g;
    }

    public final String getOStringForTest(Mem mem,int o) {
        StringBuffer r=new StringBuffer();
        int hasSolution=0;
        for(int i=0;i<x.length;i++){
            String name=mem.getName(x[i]);
            if(name.startsWith("P")){
                hasSolution++;
            }
            r.append(name);
            r.append("(");
            r.append(x[i]);
            r.append(")");
            r.append(",");
        }
        r.append(g);
        r.append("\n");
        if(hasSolution>=o){
            System.out.print(r.toString());
        }
        return r.toString();
    }
}
