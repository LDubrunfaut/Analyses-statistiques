#####################################################################
#           Script Projet Long : Analyse de données RNASeq
#####################################################################
# DUBRUNFAUT Lucie
# M2BI
# Année 2016 - 2017
#####################################################################
# Vous devez disposer dans le répertoire des fichiers suivants :
# > DUBRUNFAUT_Projet_Long
#   > genes.expected_counts
#   > sample.annotation.csv
#   > DUBRUNFAUT_script.R
#   > Ensembl_toHGNC
#       > dt1 à 5 au format txt
#   > OR
#       > OR_f1 à 14 + 52, 53, 55, 56 au format txt
#####################################################################
# Chargement des données :
#####################################################################
countData=read.table("genes.expected_counts",h=T)

# Chargement des matrices de correspondance Ensembl - HGNC
dt_1_to_10k=read.table("./Ensembl_to_HGNC/dt1.txt",sep="\t",h=T)
dt_10k_to_20k=read.table("./Ensembl_to_HGNC/dt2.txt",sep="\t",h=T)
dt_20k_to_30k=read.table("./Ensembl_to_HGNC/dt3.txt",sep="\t",h=T)
dt_30k_to_40k=read.table("./Ensembl_to_HGNC/dt4.txt",sep="\t",h=T)
dt_end=read.table("./Ensembl_to_HGNC/dt5.txt",sep="\t",h=T)

temp1=rbind(dt_1_to_10k,dt_10k_to_20k)
temp2=rbind(temp1,dt_20k_to_30k)
temp3=rbind(temp2,dt_30k_to_40k)
temp4=rbind(temp3,dt_end)
ensembl_to_hgnc=temp4[,1:2]
colnames(ensembl_to_hgnc)=c("gene_id","HGNC.ID")

# Association des données avec leurs noms HGNC
dt.all=merge(ensembl_to_hgnc,countData,by=c("gene_id"))

#####################################################################
# Mise en place des listes des gènes olfactifs :
#####################################################################
####
# Chargement de la liste des gènes olfactifs selon HGNC : 
OR1=read.table("./OR/OR_f1.txt",sep="\t",h=T)
OR2=read.table("./OR/OR_f2.txt",sep="\t",h=T)
OR3=read.table("./OR/OR_f3.txt",sep="\t",h=T)
OR4=read.table("./OR/OR_f4.txt",sep="\t",h=T)
OR5=read.table("./OR/OR_f5.txt",sep="\t",h=T)
OR6=read.table("./OR/OR_f6.txt",sep="\t",h=T)
OR7=read.table("./OR/OR_f7.txt",sep="\t",h=T)
OR8=read.table("./OR/OR_f8.txt",sep="\t",h=T)
OR9=read.table("./OR/OR_f9.txt",sep="\t",h=T)
OR10=read.table("./OR/OR_f10.txt",sep="\t",h=T)
OR11=read.table("./OR/OR_f11.txt",sep="\t",h=T)
OR12=read.table("./OR/OR_f12.txt",sep="\t",h=T)
OR13=read.table("./OR/OR_f13.txt",sep="\t",h=T)
OR14=read.table("./OR/OR_f14.txt",sep="\t",h=T)
OR51=read.table("./OR/OR_f51.txt",sep="\t",h=T)
OR52=read.table("./OR/OR_f52.txt",sep="\t",h=T)
OR55=read.table("./OR/OR_f55.txt",sep="\t",h=T)
OR56=read.table("./OR/OR_f56.txt",sep="\t",h=T)

ORs=rbind(OR1,OR2,OR3,OR4,OR5,OR6,OR7,OR8,OR9,OR10,OR11,OR12,OR13,OR14,OR51,OR52,OR55,OR56)
# 861 ORs idenfitifiés

# Sélection des gènes dans nos données appartenant à cette liste
our_OR=merge(ORs,dt.all,by=c("HGNC.ID"))
our_OR=cbind(our_OR[,1],our_OR[,12:dim(our_OR)[2]])
colnames(our_OR)[1]="HGNC.ID"
# Dans nos données il y 491 ORs identifiés

####
# Chargement de la liste de gènes mis en évidence dans la publication
# "The human olfactif transcriptome"

# On passe de HGNC à Symbol et de Symbol à ENSEMBL
symb=read.table("symbol_to_HGNC.txt",sep="\t",h=T)
symb=symb[,1:2]
ens=read.table("ensembl_to_HGNC.txt",sep="\t",h=T)
ens=ens[,1:2]
our_symb=merge(symb,dt.all,by=c("HGNC.ID"))
our_symb2=cbind(our_symb[,1],our_symb[,3:dim(our_symb)[2]])
colnames(our_symb2)[1]="HGNC.ID"
our_ens=merge(ens,dt.all,by=c("HGNC.ID"))
our_ens2=cbind(our_ens[,1],our_ens[,3:dim(our_ens)[2]])
colnames(our_ens2)[1]="HGNC.ID"

# Nouveaux gènes olfactifs
our_new_genes_olf=rbind(our_symb2,our_ens2)
# Dans nos données il y 90 des nouveaux gènes olfactifs associés à la
# publication.

# On enregistre ces deux listes :
write.table(our_OR, "ORs_in_our_data.txt", sep="\t", row.names=FALSE, na="")
write.table(our_new_genes_olf "ORs_from_publi_in_our_data.txt", sep="\t", row.names=FALSE, na="")
#####################################################################
#####################################################################
# Analyse des données :
#   --> Utilisation du package EdgeR
#####################################################################
#####################################################################
library("edgeR")

# Chargement des listes et des annotation des échantillons
liste_ORs=read.table("ORs_in_our_data.txt",sep="\t",h=T)
liste_new_genes_olf=read.table("ORs_from_publi_in_our_data.txt", sep="\t",h=T)
colData = read.csv("sample_annotation.csv", row.names=1)

# Etant donné que les identifiants HGNC ne sont pas uniques, 
# on garde les ENSEMBL comme clé unique des lignes

temp1= liste_ORs
liste_ORs = liste_ORs[,-1]
liste_ORs = liste_ORs[,-1]
rownames(liste_ORs) = temp1[,2]
colnames(liste_ORs) = colData[,1]

temp1= liste_new_genes_olf
liste_new_genes_olf = liste_new_genes_olf[,-1]
liste_new_genes_olf = liste_new_genes_olf[,-1]
rownames(liste_new_genes_olf) = temp1[,2]
colnames(liste_new_genes_olf) = colData[,1]

# Construction de l'objet des données pour EdgeR :

######################################
#### Attention vous devez ici lancer l'analyse pour chacune des deux
# listes manuelle en changeant la valeur de counts avec l'une ou 
# l'autre des variables 
######################################

dds = DGEList(counts=liste_ORs,group=factor(colData[,2]))
# dds = DGEList(counts=liste_new_genes_olf,group=factor(colData[,2]))

#####################################################################
# Filtrage des données pour ne conserver que les gènes qui varient
# dans nos données (non supervisé)
#####################################################################
# On ne conserve que les gènes ayant au moins 1 cpm dans au moins
# trois échantillons
countsPerMillion = cpm(dds)
countCheck = countsPerMillion > 1
keep = which(rowSums(countCheck) >=3)
cds = dds[keep,]

#####################################################################
# Normalisation
#####################################################################
# Les données sont déjà normalisées 
#####################################################################
# Exploration non supervisée des données :
#   --> Multidimensional Scaling
#####################################################################
dt_MDS=plotMDS(cds,method="bcv", col=as.numeric(dds$samples$group))

png(filename="MDS_liste_OR.png")
# png(filename="MDS_liste_new_genes_olf.png")
plotMDS(dt_MDS,col=as.numeric(dds$samples$group), dim.plot=c(1,2))
legend("topleft",as.character(unique(dds$samples$group)),col=1:6,pch=20)
dev.off()

#####################################################################
# Estimation de la dispersion des données
#####################################################################

d1=estimateCommonDisp(cds,verbose=T)
# Disp = 0.09518 , BCV = 0.3085

# prior.n recommandé est : 50/(nbr_d'échantillons - nbr_de_groupes)
d1=estimateTagwiseDisp(d1,prior.df=50/(24-6))

d1=estimateTrendedDisp(d1)

png(filename="Dispersion_plot_liste_OR_qCML")
# png(filename="Dispersion_plot_liste_new_genes_olf_qCML")
plotBCV(d1)
dev.off()

#####################################################################
# Analyse de l'expression différentielle
#####################################################################
####
# Valeur de cutoff pour retenir les gènes significativement exprimés
# différentiellement :
cutoff=0.20

# Analyse de l'expression différentielle en utilisant d1 :

# On souhaite comparer tous les facteurs de notre étude deux à deux :
# Remarque on compare A et B mais aussi B et A
facteurs=levels(cds$samples$group)
tTag.list=data.frame("Comparaison"="t","Genes_DE"="e")
tTag.list=tTag.list[-1,]
for(fac in facteurs){
    for(compared in facteurs){
        if(fac!=compared){
            pair = paste(paste(fac,"vs",sep=""),compared,sep="")
            topt=paste("d1.",pair,sep="")
            temp1=exactTest(d1, pair=c(fac,compared))
            # On créé un objet avec les topTags au passage
            tTag=topTags(temp1,n=nrow(temp1$table))$table
            gDE=rownames(tTag)[tTag$FDR<=cutoff]
            tTag.list=rbind(tTag.list,data.frame(Comparaison=pair,Genes_DE=paste(gDE,collapse=";")))
            # Pour chaque gène sélectionné on cherche son expression
            # différentielle :
            DE_gene=tTag[which(rownames(tTag)%in% gDE),]
            if(nrow(DE_gene)!=0){
                print(head(DE_gene))
            }
            assign(topt,DE_gene)        
        }
    }
}

# Identifiants HGNC correspondants :
tTag.list$HGNC="HGNC"
for(i in 1:nrow(tTag.list)){
    genes_id = unlist(strsplit(as.character(tTag.list[i,2]),split=";"))
    temp=subset(ensembl_to_hgnc, gene_id %in% genes_id)
    if(nrow(temp)!=0){
        print(temp)
    }
    temp1=as.vector(temp[,2])
    tTag.list[i,3]=paste(temp1,collapse=";")
}

write.table(tTag.list, "Cutoff=0,20.txt", sep="\t", row.names=FALSE, na="")






















