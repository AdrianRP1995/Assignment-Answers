

require 'net/http'
require 'json'
#We require the interacting proteins Class


#=====================================================================================================================================================


class Protein

  attr_accessor :uniprotid # UniProt ID
  attr_accessor :intactid # IntAct ID, to assess the interactions with other proteins
  attr_accessor :networkid # Network ID, to load the networks each protein is part of

  @@prot_obj = Hash.new # Class variable that will store all the data we introduce into the class
  @@prot_int_obj = Hash.new # Class variable that will store all the data we introduce into the class with an intact ID as a key


#=====================================================================================================================================================


  def initialize (params = {})

    @uniprotid = params.fetch(:uniprotid, "XXXXXX")
    @intactid = params.fetch(:intactid, nil)
    @networkid = params.fetch(:networkid, nil)

    @@prot_obj[uniprotid] = self

    if intactid
      @@prot_int_obj[intactid] = self #We only introduce the data in the hash if the protein has an Intact ID, and will use it as key
    end

  end


#=====================================================================================================================================================


  def Protein.all
     #This class method will return all the instances of the class

    return @@prot_obj

  end


#=====================================================================================================================================================


  def Protein.all_intact
    #This class method will return all the instances of the class which contains an Intact ID

    return @@prot_int_obj

  end


#=====================================================================================================================================================


  def Protein.exists(intactid)
    # This class method will check whether a protein entrance with intactid was already created

    if @@prot_int_obj.has_key?(intactid)
      return true

    else
      return false

    end

  end
#-----------------------------------------------------


#-----------------------------------------------------
  def Protein.get_intactid(geneid)
    # This class methon will look for the IntAct ID of a protein and will retieve it if it is found

    ad = URI("http://togows.org/entry/ebi-uniprot/#{geneid}/dr.json")
    res = Net::HTTP.get_response(ad)

    data = JSON.parse(res.body)

    if data[0]['IntAct']
      return data[0]['IntAct'][0][0] #If there isn't an Intact ID, the result will be nil
    else
      return nil
    end


  end


#=====================================================================================================================================================


  def Protein.create_prot (uniprotid, level, geneid = nil, intactid = nil)
    # This methond creates a protein entry, using its uniprotid and the level of interaction we are in.
    # We can introduce geneid or intact id, if we are in level 1 of interaction


    if not intactid  #We will look for the intact id if it is not provided
      intactid = Protein.get_intactid(geneid) #We call this class method, which will look for the Intact ID and will retrieve it if there is one
    end

    Protein.new(
            :uniprotid => uniprotid,
            :intactid => intactid, #If there is not an Intact ID, this parameter will be empty
            :network => nil # We first set all networks relationships to empty
            )


    level += 1 # We deep one level

  end


#=====================================================================================================================================================


end
