/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.serializer;

import com.google.inject.Inject;
import java.util.Map;
import java.util.Set;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EPackage;
import org.eclipse.xtext.Action;
import org.eclipse.xtext.Parameter;
import org.eclipse.xtext.ParserRule;
import org.eclipse.xtext.serializer.ISerializationContext;
import org.eclipse.xtext.serializer.sequencer.AbstractDelegatingSemanticSequencer;

@SuppressWarnings("all")
public abstract class AbstractElkGraphJsonSemanticSequencer extends AbstractDelegatingSemanticSequencer {

	@Inject
	private ElkGraphJsonGrammarAccess grammarAccess;
	
	@Override
	public void sequence(ISerializationContext context, EObject semanticObject) {
		EPackage epackage = semanticObject.eClass().getEPackage();
		ParserRule rule = context.getParserRule();
		Action action = context.getAssignedAction();
		Set<Parameter> parameters = context.getEnabledBooleanParameters();
		if (epackage == ElkGraphPackage.eINSTANCE)
			switch (semanticObject.eClass().getClassifierID()) {
			case ElkGraphPackage.ELK_EDGE:
				sequence_ElkEdgeSources_ElkEdgeTargets_ElkGraphElementLabels_ElkGraphElementProperties_ElkId(context, (ElkEdge) semanticObject); 
				return; 
			case ElkGraphPackage.ELK_LABEL:
				sequence_ElkGraphElementLabels_ElkGraphElementProperties_ElkId_LabelElement_ShapeElement(context, (ElkLabel) semanticObject); 
				return; 
			case ElkGraphPackage.ELK_NODE:
				sequence_ElkGraphElementLabels_ElkGraphElementProperties_ElkId_ElkNode_ElkNodeChildren_ElkNodeEdges_ElkNodePorts_ShapeElement(context, (ElkNode) semanticObject); 
				return; 
			case ElkGraphPackage.ELK_PORT:
				sequence_ElkGraphElementLabels_ElkGraphElementProperties_ElkId_ShapeElement(context, (ElkPort) semanticObject); 
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
	 *         identifier=STRING? 
	 *         (sources+=[ElkConnectableShape|STRING] sources+=[ElkConnectableShape|STRING]*)? 
	 *         (targets+=[ElkConnectableShape|STRING] targets+=[ElkConnectableShape|STRING]*)? 
	 *         (labels+=ElkLabel labels+=ElkLabel*)? 
	 *         (properties+=Property properties+=Property*)?
	 *     )+
	 */
	protected void sequence_ElkEdgeSources_ElkEdgeTargets_ElkGraphElementLabels_ElkGraphElementProperties_ElkId(ISerializationContext context, ElkEdge semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkNode returns ElkNode
	 *
	 * Constraint:
	 *     (
	 *         (identifier=STRING | x=Number | y=Number | width=Number | height=Number)? 
	 *         (children+=ElkNode children+=ElkNode*)? 
	 *         (ports+=ElkPort ports+=ElkPort*)? 
	 *         (labels+=ElkLabel labels+=ElkLabel*)? 
	 *         (containedEdges+=ElkEdge containedEdges+=ElkEdge*)? 
	 *         (properties+=Property properties+=Property*)?
	 *     )+
	 */
	protected void sequence_ElkGraphElementLabels_ElkGraphElementProperties_ElkId_ElkNode_ElkNodeChildren_ElkNodeEdges_ElkNodePorts_ShapeElement(ISerializationContext context, ElkNode semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkLabel returns ElkLabel
	 *
	 * Constraint:
	 *     (
	 *         (
	 *             text=STRING | 
	 *             identifier=STRING | 
	 *             x=Number | 
	 *             y=Number | 
	 *             width=Number | 
	 *             height=Number
	 *         )? 
	 *         (labels+=ElkLabel labels+=ElkLabel*)? 
	 *         (properties+=Property properties+=Property*)?
	 *     )+
	 */
	protected void sequence_ElkGraphElementLabels_ElkGraphElementProperties_ElkId_LabelElement_ShapeElement(ISerializationContext context, ElkLabel semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     ElkPort returns ElkPort
	 *
	 * Constraint:
	 *     (
	 *         (identifier=STRING | x=Number | y=Number | width=Number | height=Number)? 
	 *         (labels+=ElkLabel labels+=ElkLabel*)? 
	 *         (properties+=Property properties+=Property*)?
	 *     )+
	 */
	protected void sequence_ElkGraphElementLabels_ElkGraphElementProperties_ElkId_ShapeElement(ISerializationContext context, ElkPort semanticObject) {
		genericSequencer.createSequence(context, semanticObject);
	}
	
	
	/**
	 * Contexts:
	 *     Property returns ElkPropertyToValueMapEntry
	 *
	 * Constraint:
	 *     (key=PropertyKey (value=StringValue | value=NumberValue | value=BooleanValue))
	 */
	protected void sequence_Property(ISerializationContext context, Map.Entry semanticObject) {
		genericSequencer.createSequence(context, (EObject) semanticObject);
	}
	
	
}
