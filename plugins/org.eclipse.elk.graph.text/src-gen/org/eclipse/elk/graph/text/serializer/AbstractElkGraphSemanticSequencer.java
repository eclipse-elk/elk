/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.serializer;

import com.google.inject.Inject;
import java.util.Map;
import java.util.Set;
import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess;
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
public abstract class AbstractElkGraphSemanticSequencer extends AbstractDelegatingSemanticSequencer {

	@Inject
	private ElkGraphGrammarAccess grammarAccess;
	
	@Override
	public void sequence(ISerializationContext context, EObject semanticObject) {
		EPackage epackage = semanticObject.eClass().getEPackage();
		ParserRule rule = context.getParserRule();
		Action action = context.getAssignedAction();
		Set<Parameter> parameters = context.getEnabledBooleanParameters();
		if (epackage == ElkGraphPackage.eINSTANCE)
			switch (semanticObject.eClass().getClassifierID()) {
			case ElkGraphPackage.ELK_BEND_POINT:
				sequence_ElkBendPoint(context, (ElkBendPoint) semanticObject); 
				return; 
			case ElkGraphPackage.ELK_EDGE:
				sequence_EdgeLayout_ElkEdge(context, (ElkEdge) semanticObject); 
				return; 
			case ElkGraphPackage.ELK_EDGE_SECTION:
				if (rule == grammarAccess.getElkEdgeSectionRule()) {
					sequence_ElkEdgeSection(context, (ElkEdgeSection) semanticObject); 
					return; 
				}
				else if (rule == grammarAccess.getElkSingleEdgeSectionRule()) {
					sequence_ElkSingleEdgeSection(context, (ElkEdgeSection) semanticObject); 
					return; 
				}
				else break;
			case ElkGraphPackage.ELK_LABEL:
				sequence_ElkLabel_ShapeLayout(context, (ElkLabel) semanticObject); 
				return; 
			case ElkGraphPackage.ELK_NODE:
				if (rule == grammarAccess.getElkNodeRule()) {
					sequence_ElkNode_ShapeLayout(context, (ElkNode) semanticObject); 
					return; 
				}
				else if (rule == grammarAccess.getRootNodeRule()) {
					sequence_RootNode_ShapeLayout(context, (ElkNode) semanticObject); 
					return; 
				}
				else break;
			case ElkGraphPackage.ELK_PORT:
				sequence_ElkPort_ShapeLayout(context, (ElkPort) semanticObject); 
				return; 
			case ElkGraphPackage.ELK_PROPERTY_TO_VALUE_MAP_ENTRY:
				sequence_Property(context, (Map.Entry) semanticObject); 
				return; 
			}
		if (errorAcceptor != null)
			errorAcceptor.accept(diagnosticProvider.createInvalidContextOrTypeDiagnostic(semanticObject, context));
	}
	
	/**
	 * Contexts:
	 *     ElkEdge returns ElkEdge
	 *
	 * Constraint:
	 *     (
	 *         identifier=ID? 
	 *         sources+=[ElkConnectableShape|QualifiedId] 
	 *         sources+=[ElkConnectableShape|QualifiedId]* 
	 *         targets+=[ElkConnectableShape|QualifiedId] 
	 *         targets+=[ElkConnectableShape|QualifiedId]* 
	 *         (sections+=ElkSingleEdgeSection | sections+=ElkEdgeSection+)? 
	 *         properties+=Property* 
	 *         labels+=ElkLabel*
	 *     )
	 */
	protected void sequence_EdgeLayout_ElkEdge(ISerializationContext context, ElkEdge semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkBendPoint returns ElkBendPoint
	 *
	 * Constraint:
	 *     (x=Number y=Number)
	 */
	protected void sequence_ElkBendPoint(ISerializationContext context, ElkBendPoint semanticObject) {
		if (errorAcceptor != null) {
			if (transientValues.isValueTransient(semanticObject, ElkGraphPackage.Literals.ELK_BEND_POINT__X) == ValueTransient.YES)
				errorAcceptor.accept(diagnosticProvider.createFeatureValueMissing(semanticObject, ElkGraphPackage.Literals.ELK_BEND_POINT__X));
			if (transientValues.isValueTransient(semanticObject, ElkGraphPackage.Literals.ELK_BEND_POINT__Y) == ValueTransient.YES)
				errorAcceptor.accept(diagnosticProvider.createFeatureValueMissing(semanticObject, ElkGraphPackage.Literals.ELK_BEND_POINT__Y));
		}
		SequenceFeeder feeder = createSequencerFeeder(context, semanticObject);
		feeder.accept(grammarAccess.getElkBendPointAccess().getXNumberParserRuleCall_0_0(), semanticObject.getX());
		feeder.accept(grammarAccess.getElkBendPointAccess().getYNumberParserRuleCall_2_0(), semanticObject.getY());
		feeder.finish();
	}
	
	
	/**
	 * Contexts:
	 *     ElkEdgeSection returns ElkEdgeSection
	 *
	 * Constraint:
	 *     (
	 *         identifier=ID 
	 *         (outgoingSections+=[ElkEdgeSection|ID] outgoingSections+=[ElkEdgeSection|ID]*)? 
	 *         (
	 *             (incomingShape=[ElkConnectableShape|QualifiedId] | outgoingShape=[ElkConnectableShape|QualifiedId])? 
	 *             (startX=Number startY=Number)? 
	 *             (endX=Number endY=Number)?
	 *         )+ 
	 *         (bendPoints+=ElkBendPoint bendPoints+=ElkBendPoint*)? 
	 *         properties+=Property*
	 *     )
	 */
	protected void sequence_ElkEdgeSection(ISerializationContext context, ElkEdgeSection semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkLabel returns ElkLabel
	 *
	 * Constraint:
	 *     (identifier=ID? text=STRING ((x=Number y=Number) | (width=Number height=Number))* properties+=Property* labels+=ElkLabel*)
	 */
	protected void sequence_ElkLabel_ShapeLayout(ISerializationContext context, ElkLabel semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkNode returns ElkNode
	 *
	 * Constraint:
	 *     (
	 *         identifier=ID 
	 *         ((x=Number y=Number) | (width=Number height=Number))* 
	 *         properties+=Property* 
	 *         (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*
	 *     )
	 */
	protected void sequence_ElkNode_ShapeLayout(ISerializationContext context, ElkNode semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkPort returns ElkPort
	 *
	 * Constraint:
	 *     (identifier=ID ((x=Number y=Number) | (width=Number height=Number))* properties+=Property* labels+=ElkLabel*)
	 */
	protected void sequence_ElkPort_ShapeLayout(ISerializationContext context, ElkPort semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkSingleEdgeSection returns ElkEdgeSection
	 *
	 * Constraint:
	 *     (
	 *         (
	 *             (incomingShape=[ElkConnectableShape|QualifiedId] | outgoingShape=[ElkConnectableShape|QualifiedId])? 
	 *             (startX=Number startY=Number)? 
	 *             (endX=Number endY=Number)?
	 *         )+ 
	 *         (bendPoints+=ElkBendPoint bendPoints+=ElkBendPoint*)? 
	 *         properties+=Property*
	 *     )
	 */
	protected void sequence_ElkSingleEdgeSection(ISerializationContext context, ElkEdgeSection semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Property returns ElkPropertyToValueMapEntry
	 *
	 * Constraint:
	 *     (key=PropertyKey (value=StringValue | value=QualifiedIdValue | value=NumberValue | value=BooleanValue)?)
	 */
	protected void sequence_Property(ISerializationContext context, Map.Entry semanticObject) {
		genericSequencer.createSequence(context, (EObject) semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     RootNode returns ElkNode
	 *
	 * Constraint:
	 *     (
	 *         identifier=ID? 
	 *         ((x=Number y=Number) | (width=Number height=Number))* 
	 *         properties+=Property* 
	 *         (labels+=ElkLabel | ports+=ElkPort | children+=ElkNode | containedEdges+=ElkEdge)*
	 *     )
	 */
	protected void sequence_RootNode_ShapeLayout(ISerializationContext context, ElkNode semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
}
