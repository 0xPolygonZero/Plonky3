use crate::air::AluAir;
use crate::air::alu::air::FieldOperation;
use crate::air::asic::Asic;
use crate::circuit_builder::gates::event::AllEvents;
use crate::verifier::recursive_traits::RecursiveStarkGenerationConfig;
pub struct RecursiveProof<
    InputProof,
    SC: RecursiveStarkGenerationConfig<InputProof, D>,
    P,
    const D: usize,
> {
    pub add_air: AluAir<1>,
    pub sub_air: AluAir<1>,
    pub add_proof: P,
    pub sub_proof: P,
    _phantom: std::marker::PhantomData<(InputProof, SC)>,
}

pub fn prove<P, InputProof, SC, const D: usize>(
    config: &SC,
    asic: Asic<SC::Val, D>,
    all_events: AllEvents<SC::Val, D>,
    base_prove: fn(
        &SC,
        &AluAir<1>,
        p3_matrix::dense::RowMajorMatrix<SC::Val>,
        &Vec<Vec<usize>>,
    ) -> P,
) -> RecursiveProof<InputProof, SC, P, D>
where
    SC: RecursiveStarkGenerationConfig<InputProof, D>,
{
    let traces = asic.generate_trace(&all_events);

    let add_air: AluAir<1> = AluAir {
        op: FieldOperation::Add,
    };
    let sub_air: AluAir<1> = AluAir {
        op: FieldOperation::Sub,
    };

    let add_proof = base_prove(config, &add_air, traces[0].clone(), &vec![]);
    let sub_proof = base_prove(config, &sub_air, traces[1].clone(), &vec![]);

    RecursiveProof {
        add_air,
        sub_air,
        add_proof,
        sub_proof,
        _phantom: std::marker::PhantomData,
    }
}
