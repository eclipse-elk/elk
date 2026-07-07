/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.serializer;

import com.google.inject.Inject;
import java.util.Set;
import org.eclipse.elk.alg.graphviz.dot.dot.Attribute;
import org.eclipse.elk.alg.graphviz.dot.dot.AttributeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.DotPackage;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget;
import org.eclipse.elk.alg.graphviz.dot.dot.Graph;
import org.eclipse.elk.alg.graphviz.dot.dot.GraphvizModel;
import org.eclipse.elk.alg.graphviz.dot.dot.Node;
import org.eclipse.elk.alg.graphviz.dot.dot.NodeStatement;
import org.eclipse.elk.alg.graphviz.dot.dot.Port;
import org.eclipse.elk.alg.graphviz.dot.dot.Subgraph;
import org.eclipse.elk.alg.graphviz.dot.services.GraphvizDotGrammarAccess;
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
public abstract class AbstractGraphvizDotSemanticSequencer extends AbstractDelegatingSemanticSequencer {

	@Inject
	private GraphvizDotGrammarAccess grammarAccess;
	
	@Override
	public void sequence(ISerializationContext context, EObject semanticObject) {
		EPackage epackage = semanticObject.eClass().getEPackage();
		ParserRule rule = context.getParserRule();
		Action action = context.getAssignedAction();
		Set<Parameter> parameters = context.getEnabledBooleanParameters();
		if (epackage == DotPackage.eINSTANCE)
			switch (semanticObject.eClass().getClassifierID()) {
			case DotPackage.ATTRIBUTE:
				if (rule == grammarAccess.getStatementRule()
						|| rule == grammarAccess.getAttributeRule()) {
					sequence_Attribute(context, (Attribute) semanticObject); 
					return; 
				}
				else if (rule == grammarAccess.getListAttributeRule()) {
					sequence_ListAttribute(context, (Attribute) semanticObject); 
					return; 
				}
				else break;
			case DotPackage.ATTRIBUTE_STATEMENT:
				sequence_AttributeStatement(context, (AttributeStatement) semanticObject); 
				return; 
			case DotPackage.EDGE_STATEMENT:
				sequence_EdgeStatement(context, (EdgeStatement) semanticObject); 
				return; 
			case DotPackage.EDGE_TARGET:
				sequence_EdgeTarget(context, (EdgeTarget) semanticObject); 
				return; 
			case DotPackage.GRAPH:
				sequence_Graph(context, (Graph) semanticObject); 
				return; 
			case DotPackage.GRAPHVIZ_MODEL:
				sequence_GraphvizModel(context, (GraphvizModel) semanticObject); 
				return; 
			case DotPackage.NODE:
				sequence_Node(context, (Node) semanticObject); 
				return; 
			case DotPackage.NODE_STATEMENT:
				sequence_NodeStatement(context, (NodeStatement) semanticObject); 
				return; 
			case DotPackage.PORT:
				sequence_Port(context, (Port) semanticObject); 
				return; 
			case DotPackage.SUBGRAPH:
				sequence_Subgraph(context, (Subgraph) semanticObject); 
				return; 
			}
		if (errorAcceptor != null)
			errorAcceptor.accept(diagnosticProvider.createInvalidContextOrTypeDiagnostic(semanticObject, context));
	}
	
	/**
	 * Contexts:
	 *     Statement returns AttributeStatement
	 *     AttributeStatement returns AttributeStatement
	 *
	 * Constraint:
	 *     (type=AttributeType (attributes+=ListAttribute attributes+=ListAttribute*)?)
	 */
	protected void sequence_AttributeStatement(ISerializationContext context, AttributeStatement semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Statement returns Attribute
	 *     Attribute returns Attribute
	 *
	 * Constraint:
	 *     (name=DotID value=DotID)
	 */
	protected void sequence_Attribute(ISerializationContext context, Attribute semanticObject) {
		if (errorAcceptor != null) {
			if (transientValues.isValueTransient(semanticObject, DotPackage.Literals.ATTRIBUTE__NAME) == ValueTransient.YES)
				errorAcceptor.accept(diagnosticProvider.createFeatureValueMissing(semanticObject, DotPackage.Literals.ATTRIBUTE__NAME));
			if (transientValues.isValueTransient(semanticObject, DotPackage.Literals.ATTRIBUTE__VALUE) == ValueTransient.YES)
				errorAcceptor.accept(diagnosticProvider.createFeatureValueMissing(semanticObject, DotPackage.Literals.ATTRIBUTE__VALUE));
		}
		SequenceFeeder feeder = createSequencerFeeder(context, semanticObject);
		feeder.accept(grammarAccess.getAttributeAccess().getNameDotIDParserRuleCall_0_0(), semanticObject.getName());
		feeder.accept(grammarAccess.getAttributeAccess().getValueDotIDParserRuleCall_2_0(), semanticObject.getValue());
		feeder.finish();
	}
	
	
	/**
	 * Contexts:
	 *     Statement returns EdgeStatement
	 *     EdgeStatement returns EdgeStatement
	 *
	 * Constraint:
	 *     (sourceNode=Node edgeTargets+=EdgeTarget+ (attributes+=ListAttribute attributes+=ListAttribute*)?)
	 */
	protected void sequence_EdgeStatement(ISerializationContext context, EdgeStatement semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     EdgeTarget returns EdgeTarget
	 *
	 * Constraint:
	 *     (operator=EdgeOperator (targetSubgraph=Subgraph | targetnode=Node))
	 */
	protected void sequence_EdgeTarget(ISerializationContext context, EdgeTarget semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Graph returns Graph
	 *
	 * Constraint:
	 *     (strict?='strict'? type=GraphType name=DotID? statements+=Statement*)
	 */
	protected void sequence_Graph(ISerializationContext context, Graph semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     GraphvizModel returns GraphvizModel
	 *
	 * Constraint:
	 *     graphs+=Graph+
	 */
	protected void sequence_GraphvizModel(ISerializationContext context, GraphvizModel semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ListAttribute returns Attribute
	 *
	 * Constraint:
	 *     (name=DotID value=DotID?)
	 */
	protected void sequence_ListAttribute(ISerializationContext context, Attribute semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Statement returns NodeStatement
	 *     NodeStatement returns NodeStatement
	 *
	 * Constraint:
	 *     (node=Node (attributes+=ListAttribute attributes+=ListAttribute*)?)
	 */
	protected void sequence_NodeStatement(ISerializationContext context, NodeStatement semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Node returns Node
	 *
	 * Constraint:
	 *     (name=DotID port=Port?)
	 */
	protected void sequence_Node(ISerializationContext context, Node semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Port returns Port
	 *
	 * Constraint:
	 *     (name=DotID compass_pt=ID?)
	 */
	protected void sequence_Port(ISerializationContext context, Port semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Statement returns Subgraph
	 *     Subgraph returns Subgraph
	 *
	 * Constraint:
	 *     (name=ID? statements+=Statement*)
	 */
	protected void sequence_Subgraph(ISerializationContext context, Subgraph semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
}
