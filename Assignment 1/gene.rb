class Gene
  #We define the parameters
  attr_accessor :gene_id
  attr_accessor :gene_name
  attr_accessor :mutant_phenotype
  attr_accessor :linked
  @@total_gene_objects = Hash.new #A hash to store all the data we introduce in the class, and retrieve it if necessary


#======================================================================================================================================


  def initialize (params = {})
    @gene_id = params.fetch(:gene_id, "AT0G00000")
    @gene_name = params.fetch(:gene_name, "no_gene")
    @mutant_phenotype = params.fetch(:mutant_phenotype, "no_phenotype")
    @linked = params.fetch(:linked, false) #This will change if the gene is linked
#Default params, in case we don't introduce any or they aren't correct
    @@total_gene_objects[gene_id] = self #This way, we store the data to retrieve it later
  end


#======================================================================================================================================


  def Gene.all
    return @@total_gene_objects #This function calls ALL the contents of the class (Use only in cases of small amounts of data)
  end


#======================================================================================================================================


  def Gene.get_gene (id)
    if @@total_gene_objects.has_key? (id)
      return @@total_gene_objects[id]
    else
      return abort "Error: No gene with that id"
    end
    #This function returns the instance of a specific gene.
  end


#======================================================================================================================================


  def Gene.load_file (file)
    fg = File.open(file, "r")
    fg.readline
    fg.each_line do |line| #.eac_line reads every line of an object
      id,name,phenotype = line.split("\t")
      Gene.new(
        :gene_id => id,
        :gene_name  => name,
        :mutant_phenotype => phenotype
      )
    end
    fg.close
  end


#======================================================================================================================================

end
