require './seedstock.rb' #We will need this class when we want to link gene objects



class Hybridcross

  attr_accessor :parent1
  attr_accessor :parent2
  attr_accessor :f2_wild
  attr_accessor :f2_p1
  attr_accessor :f2_p2
  attr_accessor :f2_p1p2

  @@cross_obj = Array.new #This time, we have to create an array instead of a hash, since we don't have one key


#======================================================================================================================================


  def initialize (params = {})

    @parent1 = params.fetch(:parent1, "no_parent1")
    @parent2 = params.fetch(:parent2, "no_parent2") #We will use both parent objects to link both cross and seedstock classes
    @f2_wild = params.fetch(:f2_wild, "0")
    @f2_p1 = params.fetch(:f2_p1, "0")
    @f2_p2 = params.fetch(:f2_p2, "0")
    @f2_p1p2 = params.fetch(:f2_p1p2, "0")

    @@cross_obj << self #We write the data we introduce into the class in the array

  end


#======================================================================================================================================


  def Hybridcross.all
    return @@cross_obj
  end
#This function returns all the class


#======================================================================================================================================


  def Hybridcross.load_file (file)
  #This function loads the .tsv file into our class
    fc = File.open(file, "r")
    fc.readline
    fc.each_line do |line| #.each_line reads every line of an object
      parent1,parent2,f2_wild, f2p1, f2p2, f2p1p2 = line.split("\t")
      Hybridcross.new(
      :parent1 => Seedstock.get_seed_stock(parent1),
      :parent2 => Seedstock.get_seed_stock(parent2), #We link each parent to a seedstock object, which is itself linked to a gene object
      :f2_wild => f2_wild.to_i,
      :f2_p1 => f2p1.to_i,
      :f2_p2 => f2p2.to_i,
      :f2_p1p2 => f2p1p2.to_i
      )
    end
    fc.close
  end


#======================================================================================================================================


  def Hybridcross.csqlink (cross)
    # This function will calculate the chi square for each cross, and will link it if the chi square is above a said limit
    total = cross.f2_wild + cross.f2_p1 + cross.f2_p2 + cross.f2_p1p2
    e_wild = ((9*total) / 16).to_f
    e_p1 = ((3*total) / 16).to_f
    e_p2 = ((3*total) / 16).to_f
    e_p1p2 = ((1*total) / 16).to_f

    chisquare = (((cross.f2_wild - e_wild)**2) / e_wild + ((cross.f2_p1 - e_p1)**2) / e_p1 + ((cross.f2_p2 - e_p2)**2) / e_p2 + ((cross.f2_p1p2 - e_p1p2)**2) / e_p1p2).to_f
    if chisquare >= 7.82 #Maximum value of chi square. Above this number, genes will be linked. Value taken for 3 degrees of freedom, p value = 0.05

      puts "\nRecording: \n#{cross.parent1.mutant_gene_id.gene_name} is genetically linked to #{cross.parent2.mutant_gene_id.gene_name} with chisquare score #{chisquare}"
      #We call each gene name this way
      cross.parent1.mutant_gene_id.linked = cross.parent2.mutant_gene_id
      cross.parent2.mutant_gene_id.linked = cross.parent1.mutant_gene_id
      #We change the linked object from the gene class to the entry each gene is linked to

    end

  end


  #======================================================================================================================================


end
