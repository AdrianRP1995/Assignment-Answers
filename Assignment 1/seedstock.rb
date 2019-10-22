require './gene.rb' #We will need the gene class to link it to this class

class Seedstock
  attr_accessor :seedstock_id
  attr_accessor :mutant_gene_id
  attr_accessor :last_planted
  attr_accessor :storage
  attr_accessor :grams_remaining

  @@ss_objects = Hash.new #A hash to store all the data we introduce in the array


#======================================================================================================================================


  def initialize (params = {})
    @seedstock_id = params.fetch(:seedstock_id, "xxx")
    @mutant_gene_id = params.fetch(:mutant_gene_id, "AT0G00000") #We will use this variable to link both seedstock and gene classes
    @last_planted = params.fetch(:last_planted, "no_record")
    @storage = params.fetch(:storage, "no_storage")
    @grams_remaining = params.fetch(:grams_remaining, "0")

    @@ss_objects[seedstock_id] = self #This way, each time we initialize a new object (entry), we record ir to retrieve ir later

  end


#======================================================================================================================================


  def Seedstock.all
    return @@ss_objects #This function calls ALL the contents of the class
  end


#======================================================================================================================================


  def Seedstock.get_seed_stock (id)
    if @@ss_objects.has_key? (id)
      return @@ss_objects[id]
    else
      return abort "Error: No seedstock entry with that id"
    end
    #This function returns the instance of a specific seed stock entry
  end


#======================================================================================================================================


  def Seedstock.load_file (file)
    fg = File.open(file, "r")
    fg.readline
    fg.each_line do |line| #.each_line reads every line of an object
      id,mutant,last_planted, storage, grams_remaining = line.split("\t")
      Seedstock.new(
                  :seedstock_id => id,
                  :mutant_gene_id => Gene.get_gene(mutant), #We call an object of the gene class to define mutant_gene_id
                  :last_planted => last_planted,
                  :storage => storage,
                  :grams_remaining => grams_remaining.to_i #We convert the string to integer to do mathematical operations with it
                  )
    end
    fg.close
  end


#======================================================================================================================================


  def plant (amount)
    amount = amount.to_i
    if amount < @grams_remaining
      @grams_remaining -= amount
    else
      @grams_remaining = 0
      puts "\nWARNING: we have run out of Seed Stock #{@seedstock_id}" #This function substracts the amount we determine, and give us a warining if we run out of stock

    end
  end


#======================================================================================================================================


  def Seedstock.overwrite (newfile)

    if File.exists? (newfile)
      File.delete (newfile)
    end

    f = File.open(newfile, "a+")
    f.puts "Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining" #We print the header of the new file

    date = Time.now.strftime("%d/%m/%Y") #We set the time to right now
    @@ss_objects.each do |ssid, ssobject|
      ssobject.plant(7) # We call the plant method, with each class object
      ssobject.last_planted = date # We update the date

      f.puts "#{ssid}\t#{ssobject.mutant_gene_id.gene_id}\t#{date}\t#{ssobject.storage}\t#{ssobject.grams_remaining}"
      #We call our mutant_gene_id object, and then the gene_id object "inside" of it
    end


  end


#======================================================================================================================================



end
