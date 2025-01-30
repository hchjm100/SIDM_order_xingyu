void check(){

    FILE *fc =fopen("points1.txt","r");
    char line[256];
    double x;
    TH1F* h1 = new TH1F("h","h",26,-1.4,1.4);

    for(int i=0; !feof(fc) && fgets(line,256,fc);) {
       nread = sscanf(line,"%lf ", &x);
       h1->Fill(x);
    }
    fclose(fc);
    h1->Draw();
    
}
