from opentrons import protocol_api

metadata = {
    'apiLevel': '2.12',
    'protocolName': 'OCTOPUS 3.0 OT-2 Protocol: SSH-compatible, Offset Calibrated',
    'author': 'Octant',
    'description': '3.0'
}

def run(protocol: protocol_api.ProtocolContext):
    """
    OT-2 Protocol for automating the Library Pooling and Library Pool Purification 
    stages of the ExpressPlex Library Preparation kit guide.
    https://5527966.app.netsuite.com/core/media/media.nl?id=242291&c=5527966&h=OBqIoVnIJtTd8LUnt_rQUEUHw-385wvvAxRUzpb0p3QAOE7-&_xt=.pdf

    After the one-hour incubation of the ExpressPlex Reaction, the samples can 
    thereafter be processed by the OT-2 for the rest of the ExpressPlex procotol 
    up to final validation/quantification for sequencer loading.

    This specific layout allows for running on the OT-2 directly through SSH, 
    includes offset-calibration templates, and runs through all the above steps. 

    User will be prompted to type in number of plate indexes (1-8).

    Runtime is around 30 minutes.

    The following equipment/labware are assumed:
    1x  8-Channel P300 (GEN2)
    1x  8-Channel P20 (GEN2)
    1x  Magnetic Module
    1x  Opentrons 96 Filter Tip Rack 200 uL
    1x  Opentrons 96 Filter Tip Rack 20 uL
    1x  Bio-Rad 96 Well Plate 200 uL PCR
    1-2x  Applied Biosystems MicroAmp 384 Well Plate 40 uL (custom labware)
    1x  NEST 96 Deepwell Plate 2mL
    """

    ### SPECIFY True if running protocol on OT-2 command line, otherwise False if over app
    ssh_mode = True
    
    # the number of plate indexes to process: integer between 1 and 8 inclusive
    number_of_indexes = int(input('Enter number of plate indexes 1 and 8 inclusive: '))
    if number_of_indexes < 1 or number_of_indexes > 8:
        raise Exception("Error: Number of plate indexes must be between 1 and 8 inclusive.")

    protocol.set_rail_lights(True)
    print(f'** Pooling {str(number_of_indexes)} plates **')
    Run = OCTOPUS(protocol, ssh_mode)
    Run.Pool(number_of_indexes)
    Run.Purify()
    protocol.set_rail_lights(False)


class OCTOPUS:
    """The OCTOPUS class serves as a wrapper for the pooling and pool purification stages of ExpressPlex.

    Each OT-2 robot will have its own offset-calibration numbers. The ones here are only for the OT-2 used at Octant. 
    To determine the offset calibrations of one's own OT-2, Opentrons has the following guide:
    https://support.opentrons.com/en/articles/5901728-how-to-prepare-jupyter-notebook-and-opentrons_execute-for-5-0-0

    The state of SSH Mode is required since this mode is incompatible with the Opentrons App as of Server Version 5.0.2.
    """
    def __init__(self, protocol, ssh_mode=False):
        ### tip-rack definitions
            # opentrons 200 ul filtered tips
        tips_200 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '6')
            # opentrons 20 ul filtered tips
        tips_20 = protocol.load_labware('opentrons_96_filtertiprack_20ul', '5') 

        ### ADJUST LABWARE OFFSET VALUES if ssh (also later on for other labware)
        if ssh_mode:
            tips_200.set_offset(x=-0.60, y=1.00, z=-0.20)
            tips_20.set_offset(x=0.00, y=1.00, z=-0.10)

        self.tips_200 = [tips_200]
        self.tips_20 = [tips_20]

        ### pipette definitions
        self.p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=self.tips_200)
        self.p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=self.tips_20)

        ### load magnetic module
        self.mag_mod = protocol.load_module('magnetic module gen2', '9')

        ### labware specific to both Pool and Purify stages
        self.reagent_plate = self.mag_mod.load_labware('nest_96_wellplate_2ml_deep', label='Reagent Plate')
        self.dna_plate = protocol.load_labware('biorad_96_wellplate_200ul_pcr', '4', label='DNA Plate')

        ### pass down protocol api throughout object
        self.protocol = protocol

        ### pass down ssh-mode status throughout object
        self.ssh_mode = ssh_mode

    def Pool(self, num_indexes=0):
        """Pool Stage.

        The OT-2 takes over after the one-hour ExpressPlex reaction. 
        Each 96-well of samples is first consolidated into its own single column, then further consolidated into its own single well.

        Takes around 1.5 minutes per plate index.
        """
        NUM_INDEXES = num_indexes
        if NUM_INDEXES < 1 or NUM_INDEXES > 8:
            raise Exception("Error: Number of plate indexes must be between 1 and 8 inclusive.")

        ### labware specific to Pool stage
        self.reaction_plate_1 = self.protocol.load_labware('appliedbiosystems_384_wellplate_40ul', '1', label='Reaction Plate 1 (1-4)')
        self.reaction_plate_2 = self.protocol.load_labware('appliedbiosystems_384_wellplate_40ul', '2', label='Reaction Plate 2 (5-8)')

        ### ADJUST LABWARE OFFSET VALUES if ssh
        SSH_MODE = self.ssh_mode
        if SSH_MODE:
            self.reagent_plate.set_offset(x=0.10, y=1.60, z=0.00)
            self.dna_plate.set_offset(x=0.30, y=1.10, z=0.10)
            self.reaction_plate_1.set_offset(x=-0.10, y=0.80, z=0.60)
            self.reaction_plate_2.set_offset(x=-0.10, y=0.80, z=0.60)

        ### consolidate Reaction Plate(s) columns into a single column per 96-well samples
        self.set_flow_rate(2, 6, 3)
            # consolidating plates 1-4 to DNA Plate columns 1-4
        if NUM_INDEXES >= 1:
            self.p20.consolidate(2.5, [self.reaction_plate_1.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.dna_plate.columns('1'))
        if NUM_INDEXES >= 2:
            self.p20.consolidate(2.5, [self.reaction_plate_1.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.dna_plate.columns('2'))
        if NUM_INDEXES >= 3:
            self.p20.consolidate(2.5, [self.reaction_plate_1.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.dna_plate.columns('3'))
        if NUM_INDEXES >= 4:
            self.p20.consolidate(2.5, [self.reaction_plate_1.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.dna_plate.columns('4'))
            # consolidating plates 5-8 to DNA Plate columns 5-8
        if NUM_INDEXES >= 5:
            self.p20.consolidate(2.5, [self.reaction_plate_2.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.dna_plate.columns('5'))
        if NUM_INDEXES >= 6:
            self.p20.consolidate(2.5, [self.reaction_plate_2.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.dna_plate.columns('6'))
        if NUM_INDEXES >= 7:
            self.p20.consolidate(2.5, [self.reaction_plate_2.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.dna_plate.columns('7'))
        if NUM_INDEXES == 8:
            self.p20.consolidate(2.5, [self.reaction_plate_2.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.dna_plate.columns('8'))

        ### further consolidate columns into single wells on the deep well Reagent Plate
        self.set_flow_rate(25, 25, 25)
        for i in range(NUM_INDEXES):
            # pick up SINGLE tips with multichannel starting from bottom left of tip box
            self.p300.pick_up_tip(self.tips_200[0].rows()[7-i][0])
            # aspirate column
            for j in range(8):
                self.p300.aspirate(25, self.dna_plate.rows()[j][i].bottom(z=0.1))
            self.p300.dispense(200, self.reagent_plate.rows()[i][0])
            self.p300.drop_tip()

            
    def Purify(self):
        """Purify Stage.

        Pooled libraries are concentrated and purified.
        Ensure MAGwise beads are fully resuspended as this protocol does not mix beads beforehand.
        Supernatant waste is largely discarded into the deep-well plate by default.
        Magnets are left risen at the end in case DNA eluate needs to be manually transferred out.

        Takes around 26 minutes for any number of plates.
        """

        ### ADJUST LABWARE OFFSET VALUES if ssh
        SSH_MODE = self.ssh_mode
        if SSH_MODE:
            self.reagent_plate.set_offset(x=0.10, y=1.60, z=0.00)
            self.dna_plate.set_offset(x=0.30, y=1.10, z=0.10)

        ### Reagent Plate reagent locations
        DNA_pool = self.reagent_plate.wells_by_name()['A1']
        MAGwise_Beads = self.reagent_plate.wells_by_name()['A2']
        Ethanol = self.reagent_plate.wells_by_name()['A5']
        MAGwise_Beads_discard = self.reagent_plate.wells_by_name()['A3']
        Ethanol_discard = self.reagent_plate.wells_by_name()['A3']
        Tris_HCl = self.reagent_plate.wells_by_name()['A6']

        ### add and incubate MAGwise Beads with pooled samples
        self.p300.pick_up_tip()
            # adding beads to pool
        self.set_flow_rate(50, 50, 50)
        self.p300.aspirate(150, MAGwise_Beads)
        self.p300.dispense(150, DNA_pool)
            # mix beads well with pool and incubate
        self.set_flow_rate(48*4, 48*6, 48)
        self.p300.mix(3, 200, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48*4, 48*4, 48)
        for i in range(15):
            self.p300.aspirate(200, DNA_pool.bottom(z=1.0))
            self.p300.dispense(200, DNA_pool.bottom(z=10.0))
        self.set_flow_rate(48*4, 48*6, 48/6)
        self.p300.mix(3, 200, DNA_pool.bottom(z=1.0))
        self.p300.move_to(DNA_pool.bottom(z=10.0))
        self.protocol.delay(seconds=5.0)
        self.p300.blow_out(DNA_pool.bottom(z=10.0))
        self.p300.drop_tip()
        self.protocol.delay(minutes=5)
        self.mag_mod.engage()
        self.protocol.delay(minutes=5)

        ### discard supernatant in stages to minimize bead loss
        self.set_flow_rate(48/2, 48, 48)
        self.p300.pick_up_tip()
        self.p300.aspirate(200, DNA_pool.bottom(z=0.5))
        self.p300.dispense(200, MAGwise_Beads_discard)
        self.p300.move_to(MAGwise_Beads_discard.top(z=10.0))
        self.protocol.delay(minutes=2)
        self.p300.aspirate(150, DNA_pool.bottom(z=0.5))
        self.p300.dispense(150, MAGwise_Beads_discard)
        self.p300.drop_tip()
        self.set_flow_rate(20, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(20, DNA_pool.bottom())
        self.p20.drop_tip()

        ### wash with 80% ethanol twice
        self.set_flow_rate(48, 48, 48)
        self.p300.pick_up_tip()
        for i in range(1):
            self.p300.aspirate(200, Ethanol)
            self.p300.dispense(200, DNA_pool.bottom(z=10.0))
        self.p300.move_to(DNA_pool.top(z=10.0))
        self.protocol.delay(seconds=40)
        self.discard_supernatant(200, DNA_pool, Ethanol_discard.top(z=-10.0))
            # second time
        self.set_flow_rate(48, 48, 48)
        self.p300.pick_up_tip()
        for i in range(1):
            self.p300.aspirate(200, Ethanol)
            self.p300.dispense(200, DNA_pool.bottom(z=10.0))
        self.p300.move_to(DNA_pool.top(z=10.0))
        self.protocol.delay(seconds=40)
        self.discard_supernatant(200, DNA_pool, Ethanol_discard.top(z=-10.0))

        ### discard residual ethanol
        self.set_flow_rate(5, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(10, DNA_pool.bottom(z=-0.5))
        self.p20.drop_tip()

        ### resuspend pellet well in 30 ul
        self.mag_mod.disengage()
        self.set_flow_rate(30, 30, 30)
        self.p300.pick_up_tip()
        self.p300.aspirate(30, Tris_HCl.bottom(z=0.5))
        self.p300.dispense(30, DNA_pool)
        self.set_flow_rate(48/2, 48/2, 48/2)
        self.p300.mix(3, 25, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48, 48, 48)
        for i in range(20):
            self.p300.aspirate(20, DNA_pool.bottom(z=1.0))
            self.p300.dispense(20, DNA_pool.bottom(z=4.0))
        self.set_flow_rate(48/2, 48/2, 48/2)
        self.p300.mix(3, 25, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48, 48*2, 48/2)
        self.p300.mix(10, 18, DNA_pool.bottom(z=1.0))
        self.p300.blow_out(DNA_pool.bottom(z=4.0))
        self.p300.drop_tip()
        self.protocol.delay(minutes=5)
        self.mag_mod.engage()
        self.protocol.delay(minutes=2)

        ### end of library pool purification
        self.set_flow_rate(20, 20, 5)
        self.p20.pick_up_tip()
        self.p20.aspirate(20, DNA_pool)
        self.p20.dispense(20, self.dna_plate.wells_by_name()['A12'])
        self.set_flow_rate(5, 5, 5)
        self.p20.aspirate(5, DNA_pool.bottom())
        self.p20.dispense(5, self.dna_plate.wells_by_name()['A12'])
        self.set_flow_rate(1, 5, 5)
        self.p20.blow_out(self.dna_plate.wells_by_name()['A12'].top())
        self.p20.drop_tip()
        self.protocol.home()
        if SSH_MODE:
            input('Press ENTER to lower magnets and complete protocol.')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Resume to lower magnets and complete protocol.')
        self.mag_mod.disengage()
        self.protocol.home()


    def set_flow_rate(self, aspirate=96, dispense=96, blow_out=96):
        """Set the speeds (in uL/s) for all pipettes."""
        self.p300.flow_rate.aspirate = aspirate
        self.p300.flow_rate.dispense = dispense
        self.p300.flow_rate.blow_out = blow_out
        self.p20.flow_rate.aspirate = aspirate
        self.p20.flow_rate.dispense = dispense
        self.p20.flow_rate.blow_out = blow_out

    def discard_supernatant(self, vol, well, trash=None):
        """Discard liquid from wells using the P300 multi-channel."""
        if trash is None:
            trash = self.protocol.fixed_trash['A1']
        if not self.p300.has_tip:
            self.p300.pick_up_tip()
        self.set_flow_rate(48/2, 48, 48)
        while vol > 0:
            if vol < 200:
                self.p300.aspirate(vol, well.bottom(z=0.5))
                self.p300.dispense(vol, trash)
            else:
                self.p300.aspirate(200, well.bottom(z=0.5))
                self.p300.dispense(200, trash)
            vol = vol - 200
        self.p300.drop_tip()
        self.set_flow_rate(48, 48, 48)

    def discard_residual(self, well):
        """Discard remaining liquid from wells using the P20 multi-channel."""
        self.set_flow_rate(1, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(5, well.bottom())
        self.set_flow_rate(1/2, 48, 48)
        self.p20.aspirate(2, well.bottom(z=-0.5))
        self.p20.drop_tip()
        self.set_flow_rate(48, 48, 48)
