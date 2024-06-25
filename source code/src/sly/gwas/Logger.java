package sly.gwas;

public class Logger {
    public final static int INFO=0;
    public final static int DEBUG=1;
    private static int level=INFO;
    public final static void info(String s){
        if(level>=INFO){
            System.out.print("info: ");
            System.out.print(s);
        }
    }
    public final static void debug(String s){
        if(level>=DEBUG){
            System.out.print("debug: ");
            System.out.print(s);
        }
    }
    public final static void error(String s){
        System.out.print("error: ");
        System.out.print(s);
    }
    public final static void setLevel(int l){
        level=l;
    }
}
