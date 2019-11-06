require 'net/http'
require 'json'
require './protein.rb'
#We will use a protein class in case we want to access the protein ids


class Gene

  attr_accessor :geneid
  attr_accessor :uniprotid
  attr_accessor :keggid # We will create a hash to contain KEGG annotations (key: KEGG Pathway ID, value: KEGG Pathway Name)
  attr_accessor :goid # We will create a hash to contain GO annotations (key: GO ID, value: GO Term name)

  @@gene_obj = Hash.new # Class variable that will store all the data we introduce into the class


#=====================================================================================================================================================


  def initialize (params = {})

    @geneid = params.fetch(:geneid, "AT0G00000")
    @uniprotid = params.fetch(:uniprotid, "XXXXXX")

    #We will only use these entries in case they have a KEGG or GO ID
    @keggid = params.fetch(:keggid, Hash.new) # The hash is empty by default
    @goid = params.fetch(:goid, Hash.new) # The hash is empty by default


    @@gene_obj[uniprotid] = self # When we introduce data into the class, it will be recorded here

  end


#=====================================================================================================================================================


  def Gene.get_geneid
    # This class method will return all the instances of the class

    return @@gene_obj

  end


#=====================================================================================================================================================


  def Gene.get_uniprotid(geneid)
    # This class method will use the Gene ID to get the corresponding Uniprot ID

    ad = URI("http://togows.org/entry/ebi-uniprot/#{geneid}/entry_id.json")
    res = Net::HTTP.get_response(ad)

    data = JSON.parse(res.body)

    return data[0]  #This address has one entry only, which is the Uniprot ID

  end


#=====================================================================================================================================================


  def Gene.load_file(file)
    #This class method will open the file with the list of genes, read it, and introduce the data into the Gene Class
    #It will also call the Gene.get_uniprotid method to get the Uniprot ID.

    f = File.open(file, "r")

    f.each_line do |line|

      line.delete!("\n") #We have to remove the \n in each line

      uniprotid = Gene.get_uniprotid(line)

      Gene.new(
            :geneid => line,
            :uniprotid => uniprotid
            )

      Protein.create_prot(uniprotid, 0, geneid = line) # We call this method from Protein class to create a protein object

    end

    f.close

  end


#=====================================================================================================================================================


  def annotate
    # If there is a KEGG and/or GO annotation, this function will load it into the class

    adK = URI("http://togows.org/entry/kegg-genes/ath:#{Gene.geneid}/pathways.json")
    adG = URI("http://togows.org/entry/ebi-uniprot/#{Gene.geneid}/dr.json")

    resK = Net::HTTP.get_response(adK)
    resG = Net::HTTP.get_response(adG)

    dataK = JSON.parse(resK.body)
    dataG = JSON.parse(resG.body)


    # KEGG Pathways annotation
    if dataK[0]
      dataK[0].each do |path_id, path_name|
        Gene.keggid[path_id] = path_name # We use the pathway id as a key in the hash, and the pathway name as a value
    end


    # GO annotation
    if dataG[0]["GO"]
      dataG[0]["GO"].each do |num|
        Gene.goid[num[0]] = num[1].sub(/P:/, "") # We use the GO id as a key in the hash, and the term name as a value
        end
      end
    end


  end

  #=====================================================================================================================================================


end
