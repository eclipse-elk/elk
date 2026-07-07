/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.debug.grandom.serializer;

import com.google.inject.Inject;
import java.util.Set;
import org.eclipse.elk.core.debug.grandom.gRandom.Configuration;
import org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity;
import org.eclipse.elk.core.debug.grandom.gRandom.Edges;
import org.eclipse.elk.core.debug.grandom.gRandom.Flow;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;
import org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy;
import org.eclipse.elk.core.debug.grandom.gRandom.Nodes;
import org.eclipse.elk.core.debug.grandom.gRandom.Ports;
import org.eclipse.elk.core.debug.grandom.gRandom.RandGraph;
import org.eclipse.elk.core.debug.grandom.gRandom.Size;
import org.eclipse.elk.core.debug.grandom.services.GRandomGrammarAccess;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EPackage;
import org.eclipse.xtext.Action;
import org.eclipse.xtext.Parameter;
import org.eclipse.xtext.ParserRule;
import org.eclipse.xtext.serializer.ISerializationContext;
import org.eclipse.xtext.serializer.acceptor.SequenceFeeder;
import org.eclipse.xtext.serializer.sequencer.AbstractDelegatingSemanticSequencer;
import org.eclipse.xtext.serializer.sequencer.ITransientValueService.ValueTransient;

@SuppressWarnings("all")
public class GRandomSemanticSequencer extends AbstractDelegatingSemanticSequencer {

	@Inject
	private GRandomGrammarAccess grammarAccess;
	
	@Override
	public void sequence(ISerializationContext context, EObject semanticObject) {
		EPackage epackage = semanticObject.eClass().getEPackage();
		ParserRule rule = context.getParserRule();
		Action action = context.getAssignedAction();
		Set<Parameter> parameters = context.getEnabledBooleanParameters();
		if (epackage == GRandomPackage.eINSTANCE)
			switch (semanticObject.eClass().getClassifierID()) {
			case GRandomPackage.CONFIGURATION:
				sequence_Configuration(context, (Configuration) semanticObject); 
				return; 
			case GRandomPackage.DOUBLE_QUANTITY:
				sequence_DoubleQuantity(context, (DoubleQuantity) semanticObject); 
				return; 
			case GRandomPackage.EDGES:
				sequence_Edges(context, (Edges) semanticObject); 
				return; 
			case GRandomPackage.FLOW:
				sequence_Flow(context, (Flow) semanticObject); 
				return; 
			case GRandomPackage.HIERARCHY:
				sequence_Hierarchy(context, (Hierarchy) semanticObject); 
				return; 
			case GRandomPackage.NODES:
				sequence_Nodes(context, (Nodes) semanticObject); 
				return; 
			case GRandomPackage.PORTS:
				sequence_Ports(context, (Ports) semanticObject); 
				return; 
			case GRandomPackage.RAND_GRAPH:
				sequence_RandGraph(context, (RandGraph) semanticObject); 
				return; 
			case GRandomPackage.SIZE:
				sequence_Size(context, (Size) semanticObject); 
				return; 
			}
		if (errorAcceptor != null)
			errorAcceptor.accept(diagnosticProvider.createInvalidContextOrTypeDiagnostic(semanticObject, context));
	}
	
	/**
	 * Contexts:
	 *     Configuration returns Configuration
	 *
	 * Constraint:
	 *     (
	 *         samples=INT 
	 *         form=Form 
	 *         (
	 *             (
	 *                 nodes=Nodes | 
	 *                 edges=Edges | 
	 *                 hierarchy=Hierarchy | 
	 *                 seed=Integer | 
	 *                 format=Formats | 
	 *                 filename=STRING
	 *             )? 
	 *             (mW?='maxWidth' maxWidth=Integer)? 
	 *             (mD?='maxDegree' maxDegree=Integer)? 
	 *             (pF?='partitionFraction' fraction=DoubleQuantity)?
	 *         )+
	 *     )
	 */
	protected void sequence_Configuration(ISerializationContext context, Configuration semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     DoubleQuantity returns DoubleQuantity
	 *
	 * Constraint:
	 *     (quant=Double | (min=Double minMax?='to' max=Double) | (mean=Double gaussian?=Pm stddv=Double))
	 */
	protected void sequence_DoubleQuantity(ISerializationContext context, DoubleQuantity semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Edges returns Edges
	 *
	 * Constraint:
	 *     (
	 *         (density?='density' | total?='total' | relative?='relative' | outbound?='outgoing') 
	 *         nEdges=DoubleQuantity 
	 *         (labels?='labels' | selfLoops?='self loops')*
	 *     )
	 */
	protected void sequence_Edges(ISerializationContext context, Edges semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Flow returns Flow
	 *
	 * Constraint:
	 *     (flowType=FlowType side=Side amount=DoubleQuantity)
	 */
	protected void sequence_Flow(ISerializationContext context, Flow semanticObject) {
		if (errorAcceptor != null) {
			if (transientValues.isValueTransient(semanticObject, GRandomPackage.Literals.FLOW__FLOW_TYPE) == ValueTransient.YES)
				errorAcceptor.accept(diagnosticProvider.createFeatureValueMissing(semanticObject, GRandomPackage.Literals.FLOW__FLOW_TYPE));
			if (transientValues.isValueTransient(semanticObject, GRandomPackage.Literals.FLOW__SIDE) == ValueTransient.YES)
				errorAcceptor.accept(diagnosticProvider.createFeatureValueMissing(semanticObject, GRandomPackage.Literals.FLOW__SIDE));
			if (transientValues.isValueTransient(semanticObject, GRandomPackage.Literals.FLOW__AMOUNT) == ValueTransient.YES)
				errorAcceptor.accept(diagnosticProvider.createFeatureValueMissing(semanticObject, GRandomPackage.Literals.FLOW__AMOUNT));
		}
		SequenceFeeder feeder = createSequencerFeeder(context, semanticObject);
		feeder.accept(grammarAccess.getFlowAccess().getFlowTypeFlowTypeEnumRuleCall_0_0(), semanticObject.getFlowType());
		feeder.accept(grammarAccess.getFlowAccess().getSideSideEnumRuleCall_1_0(), semanticObject.getSide());
		feeder.accept(grammarAccess.getFlowAccess().getAmountDoubleQuantityParserRuleCall_3_0(), semanticObject.getAmount());
		feeder.finish();
	}
	
	
	/**
	 * Contexts:
	 *     Hierarchy returns Hierarchy
	 *
	 * Constraint:
	 *     (levels=DoubleQuantity | edges=DoubleQuantity | numHierarchNodes=DoubleQuantity | crossHierarchRel=DoubleQuantity)*
	 */
	protected void sequence_Hierarchy(ISerializationContext context, Hierarchy semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Nodes returns Nodes
	 *
	 * Constraint:
	 *     (nNodes=DoubleQuantity (ports=Ports | labels?=Labels | size=Size | removeIsolated?='remove isolated')*)
	 */
	protected void sequence_Nodes(ISerializationContext context, Nodes semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Ports returns Ports
	 *
	 * Constraint:
	 *     (labels?=Labels | reUse=DoubleQuantity | size=Size | constraint=ConstraintType | flow+=Flow)*
	 */
	protected void sequence_Ports(ISerializationContext context, Ports semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     RandGraph returns RandGraph
	 *
	 * Constraint:
	 *     configs+=Configuration+
	 */
	protected void sequence_RandGraph(ISerializationContext context, RandGraph semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Size returns Size
	 *
	 * Constraint:
	 *     (height=DoubleQuantity | width=DoubleQuantity)*
	 */
	protected void sequence_Size(ISerializationContext context, Size semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
}
