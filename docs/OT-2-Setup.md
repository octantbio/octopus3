**OT-2 Setup**
==============

Below gives some overview of the layouts involved. Refer to the [OT-2 script itself](/ot2/octopus3-ot2-pool-purify.py) for more of the specifics, as the comments walk through the steps of the pooling and purification stage. It is recommended to run it [via SSH access to the OT-2](https://support.opentrons.com/s/article/Setting-up-SSH-access-to-your-OT-2). Otherwise, accordingly adjust the variables in the script with a text editor.

The OT-2 assists in the pooling and purification ([step 2.2](/docs/Bench-Protocol.md#22-pooling)) of 96 to 8x96 samples per run. The user is prompted at the beginning to enter the number of plate indexes to process (1 to 8). Each plate-index of samples are pooled and purified separately, with the order following this quadrant convention (A1, A2, B1, B2):

![384-well quadrant convention A1, A2, B1, B2](/img/quadrants-384-well-plate.png)

 Pooling takes about 1.5 minutes per plate index. Purification takes about 26 minutes for any number of plates. The magstand is left engaged at the end in case the user wants to manually utilize it further (a prompt will ask when to disengage).

 The Reagent Plate and 200 µL tip rack are utilized in such a way that it is possible to rotate them 180˚ for reuse, as only at most half is spent per run.

### **Required Equipment**
Refer to the [Opentrons Labware Library](https://labware.opentrons.com/?category=wellPlate) for the specific labware names.

- 1x  8-Channel P300 (GEN2)
- 1x  8-Channel P20 (GEN2)
- 1x  Magnetic Module
- 1x  Opentrons 96 Filter Tip Rack 200 µL
- 1x  Opentrons 96 Filter Tip Rack 20 µL
- 1x  Bio-Rad 96 Well Plate 200 µL PCR
- 1x  NEST 96 Deepwell Plate 2mL
- 1-2x  Applied Biosystems MicroAmp 384 Well Plate 40 µL (custom labware)

We don't use the API-provided definition (appliedbiosystemsmicroamp_384_wellplate_40ul) for the 384-well plate but instead created a custom labware definition for it that can be found [here](/ot2/appliedbiosystems_384_wellplate_40ul.json).

The 384-well plates we've been using have [catalog number A36931](https://www.thermofisher.com/order/catalog/product/A36931).

### **Final Deck Layout**

![OT-2 Deck Layout](/img/ot2-deck-layout.png)

An easy way to remember the tip-rack order is that the smaller 20 µL is placed on the smaller number deck #5 than the larger 200 µL placed on #6.

Keep deck #3 empty AT ALL TIMES! The P300 multichannel will singly pick up tips and needs full clearance in this area.

### **Reagent Plate Layout**

![OT-2 Reagent Plate](/img/ot2-reagent-plate.png)

The volumes indicated are the suggested minimum volumes to fill the Reagent Plate with for reliable transfers. The actual volumes used are about 100 µL less than what is indicated.

It is recommended to only fill the plate with reagents right before running the OT-2 to ensure MAGwise beads are resuspended (OT-2 does not mix beforehand) and prevent excess ethanol evaporation.

The other half of this plate can be used for the next run by rotating the plate 180˚.

### **DNA Plate Layout**

![OT-2 DNA Plate](/img/ot2-dna-plate.png)

Do not fill the DNA Plate with anything beforehand. Simply place a new plate on deck #4. The layout above is just to provide a visual of where the samples will eventually be on the plate.

