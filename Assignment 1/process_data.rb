require './gene.rb'
require './seedstock.rb'
require './cross.rb'

#======================================================================================================================================
#We load the data on each class:
Gene.load_file ('gene_information.tsv')

Seedstock.load_file ('seed_stock_data.tsv')

Hybridcross.load_file ('cross_data.tsv')

#======================================================================================================================================
#We simulate to plant 7 grams of each seed stock, create the new file and create a report on the empty seed stocks:
Seedstock.overwrite ('new_stock_file.tsv')

#======================================================================================================================================
#We retrieve all the hybridcross data, calculate chi square for each cross and link the genes if necessary. This methon create a report:
Hybridcross.all.each do |instance|
  Hybridcross.csqlink (instance)
end

#======================================================================================================================================
#We report al the linkages:
puts "\nFinal report:\n"

Gene.all.each do |id, gene|
  if gene.linked
    puts "#{gene.gene_name} is linked to #{gene.linked.gene_name}"
  end
end

#======================================================================================================================================
puts "\n\nGoodbye!"
