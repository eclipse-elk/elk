/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.ide.contentassist.antlr;

import com.google.common.collect.ImmutableMap;
import com.google.inject.Inject;
import com.google.inject.Singleton;
import java.util.Map;
import org.eclipse.elk.graph.json.text.ide.contentassist.antlr.internal.InternalElkGraphJsonParser;
import org.eclipse.elk.graph.json.text.services.ElkGraphJsonGrammarAccess;
import org.eclipse.xtext.AbstractElement;
import org.eclipse.xtext.ide.editor.contentassist.antlr.AbstractContentAssistParser;

public class ElkGraphJsonParser extends AbstractContentAssistParser {

	@Singleton
	public static final class NameMappings {
		
		private final Map<AbstractElement, String> mappings;
		
		@Inject
		public NameMappings(ElkGraphJsonGrammarAccess grammarAccess) {
			ImmutableMap.Builder<AbstractElement, String> builder = ImmutableMap.builder();
			init(builder, grammarAccess);
			this.mappings = builder.build();
		}
		
		public String getRuleName(AbstractElement element) {
			return mappings.get(element);
		}
		
		private static void init(ImmutableMap.Builder<AbstractElement, String> builder, ElkGraphJsonGrammarAccess grammarAccess) {
			builder.put(grammarAccess.getNodeElementAccess().getAlternatives(), "rule__NodeElement__Alternatives");
			builder.put(grammarAccess.getPortElementAccess().getAlternatives(), "rule__PortElement__Alternatives");
			builder.put(grammarAccess.getLabelElementAccess().getAlternatives(), "rule__LabelElement__Alternatives");
			builder.put(grammarAccess.getEdgeElementAccess().getAlternatives(), "rule__EdgeElement__Alternatives");
			builder.put(grammarAccess.getShapeElementAccess().getAlternatives(), "rule__ShapeElement__Alternatives");
			builder.put(grammarAccess.getPropertyAccess().getAlternatives_2(), "rule__Property__Alternatives_2");
			builder.put(grammarAccess.getPropertyKeyAccess().getAlternatives(), "rule__PropertyKey__Alternatives");
			builder.put(grammarAccess.getNumberValueAccess().getAlternatives(), "rule__NumberValue__Alternatives");
			builder.put(grammarAccess.getBooleanValueAccess().getAlternatives(), "rule__BooleanValue__Alternatives");
			builder.put(grammarAccess.getNumberAccess().getAlternatives(), "rule__Number__Alternatives");
			builder.put(grammarAccess.getJsonMemberAccess().getAlternatives_0(), "rule__JsonMember__Alternatives_0");
			builder.put(grammarAccess.getJsonValueAccess().getAlternatives(), "rule__JsonValue__Alternatives");
			builder.put(grammarAccess.getKeyChildrenAccess().getAlternatives(), "rule__KeyChildren__Alternatives");
			builder.put(grammarAccess.getKeyPortsAccess().getAlternatives(), "rule__KeyPorts__Alternatives");
			builder.put(grammarAccess.getKeyLabelsAccess().getAlternatives(), "rule__KeyLabels__Alternatives");
			builder.put(grammarAccess.getKeyEdgesAccess().getAlternatives(), "rule__KeyEdges__Alternatives");
			builder.put(grammarAccess.getKeyLayoutOptionsAccess().getAlternatives(), "rule__KeyLayoutOptions__Alternatives");
			builder.put(grammarAccess.getKeyIdAccess().getAlternatives(), "rule__KeyId__Alternatives");
			builder.put(grammarAccess.getKeyXAccess().getAlternatives(), "rule__KeyX__Alternatives");
			builder.put(grammarAccess.getKeyYAccess().getAlternatives(), "rule__KeyY__Alternatives");
			builder.put(grammarAccess.getKeyWidthAccess().getAlternatives(), "rule__KeyWidth__Alternatives");
			builder.put(grammarAccess.getKeyHeightAccess().getAlternatives(), "rule__KeyHeight__Alternatives");
			builder.put(grammarAccess.getKeySourcesAccess().getAlternatives(), "rule__KeySources__Alternatives");
			builder.put(grammarAccess.getKeyTargetsAccess().getAlternatives(), "rule__KeyTargets__Alternatives");
			builder.put(grammarAccess.getKeyTextAccess().getAlternatives(), "rule__KeyText__Alternatives");
			builder.put(grammarAccess.getElkNodeAccess().getGroup(), "rule__ElkNode__Group__0");
			builder.put(grammarAccess.getElkNodeAccess().getGroup_2(), "rule__ElkNode__Group_2__0");
			builder.put(grammarAccess.getElkNodeAccess().getGroup_2_1(), "rule__ElkNode__Group_2_1__0");
			builder.put(grammarAccess.getNodeElementAccess().getGroup_2(), "rule__NodeElement__Group_2__0");
			builder.put(grammarAccess.getNodeElementAccess().getGroup_3(), "rule__NodeElement__Group_3__0");
			builder.put(grammarAccess.getNodeElementAccess().getGroup_4(), "rule__NodeElement__Group_4__0");
			builder.put(grammarAccess.getNodeElementAccess().getGroup_5(), "rule__NodeElement__Group_5__0");
			builder.put(grammarAccess.getNodeElementAccess().getGroup_6(), "rule__NodeElement__Group_6__0");
			builder.put(grammarAccess.getElkPortAccess().getGroup(), "rule__ElkPort__Group__0");
			builder.put(grammarAccess.getElkPortAccess().getGroup_2(), "rule__ElkPort__Group_2__0");
			builder.put(grammarAccess.getPortElementAccess().getGroup_2(), "rule__PortElement__Group_2__0");
			builder.put(grammarAccess.getPortElementAccess().getGroup_3(), "rule__PortElement__Group_3__0");
			builder.put(grammarAccess.getElkLabelAccess().getGroup(), "rule__ElkLabel__Group__0");
			builder.put(grammarAccess.getElkLabelAccess().getGroup_2(), "rule__ElkLabel__Group_2__0");
			builder.put(grammarAccess.getLabelElementAccess().getGroup_2(), "rule__LabelElement__Group_2__0");
			builder.put(grammarAccess.getLabelElementAccess().getGroup_3(), "rule__LabelElement__Group_3__0");
			builder.put(grammarAccess.getLabelElementAccess().getGroup_4(), "rule__LabelElement__Group_4__0");
			builder.put(grammarAccess.getElkEdgeAccess().getGroup(), "rule__ElkEdge__Group__0");
			builder.put(grammarAccess.getElkEdgeAccess().getGroup_2(), "rule__ElkEdge__Group_2__0");
			builder.put(grammarAccess.getEdgeElementAccess().getGroup_1(), "rule__EdgeElement__Group_1__0");
			builder.put(grammarAccess.getEdgeElementAccess().getGroup_2(), "rule__EdgeElement__Group_2__0");
			builder.put(grammarAccess.getEdgeElementAccess().getGroup_3(), "rule__EdgeElement__Group_3__0");
			builder.put(grammarAccess.getEdgeElementAccess().getGroup_4(), "rule__EdgeElement__Group_4__0");
			builder.put(grammarAccess.getElkEdgeSourcesAccess().getGroup(), "rule__ElkEdgeSources__Group__0");
			builder.put(grammarAccess.getElkEdgeSourcesAccess().getGroup_1(), "rule__ElkEdgeSources__Group_1__0");
			builder.put(grammarAccess.getElkEdgeSourcesAccess().getGroup_1_1(), "rule__ElkEdgeSources__Group_1_1__0");
			builder.put(grammarAccess.getElkEdgeTargetsAccess().getGroup(), "rule__ElkEdgeTargets__Group__0");
			builder.put(grammarAccess.getElkEdgeTargetsAccess().getGroup_1(), "rule__ElkEdgeTargets__Group_1__0");
			builder.put(grammarAccess.getElkEdgeTargetsAccess().getGroup_1_1(), "rule__ElkEdgeTargets__Group_1_1__0");
			builder.put(grammarAccess.getElkIdAccess().getGroup(), "rule__ElkId__Group__0");
			builder.put(grammarAccess.getElkNodeChildrenAccess().getGroup(), "rule__ElkNodeChildren__Group__0");
			builder.put(grammarAccess.getElkNodeChildrenAccess().getGroup_1(), "rule__ElkNodeChildren__Group_1__0");
			builder.put(grammarAccess.getElkNodeChildrenAccess().getGroup_1_1(), "rule__ElkNodeChildren__Group_1_1__0");
			builder.put(grammarAccess.getElkNodePortsAccess().getGroup(), "rule__ElkNodePorts__Group__0");
			builder.put(grammarAccess.getElkNodePortsAccess().getGroup_1(), "rule__ElkNodePorts__Group_1__0");
			builder.put(grammarAccess.getElkNodePortsAccess().getGroup_1_1(), "rule__ElkNodePorts__Group_1_1__0");
			builder.put(grammarAccess.getElkNodeEdgesAccess().getGroup(), "rule__ElkNodeEdges__Group__0");
			builder.put(grammarAccess.getElkNodeEdgesAccess().getGroup_1(), "rule__ElkNodeEdges__Group_1__0");
			builder.put(grammarAccess.getElkNodeEdgesAccess().getGroup_1_1(), "rule__ElkNodeEdges__Group_1_1__0");
			builder.put(grammarAccess.getElkGraphElementLabelsAccess().getGroup(), "rule__ElkGraphElementLabels__Group__0");
			builder.put(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1(), "rule__ElkGraphElementLabels__Group_1__0");
			builder.put(grammarAccess.getElkGraphElementLabelsAccess().getGroup_1_1(), "rule__ElkGraphElementLabels__Group_1_1__0");
			builder.put(grammarAccess.getElkGraphElementPropertiesAccess().getGroup(), "rule__ElkGraphElementProperties__Group__0");
			builder.put(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1(), "rule__ElkGraphElementProperties__Group_1__0");
			builder.put(grammarAccess.getElkGraphElementPropertiesAccess().getGroup_1_1(), "rule__ElkGraphElementProperties__Group_1_1__0");
			builder.put(grammarAccess.getShapeElementAccess().getGroup_0(), "rule__ShapeElement__Group_0__0");
			builder.put(grammarAccess.getShapeElementAccess().getGroup_1(), "rule__ShapeElement__Group_1__0");
			builder.put(grammarAccess.getShapeElementAccess().getGroup_2(), "rule__ShapeElement__Group_2__0");
			builder.put(grammarAccess.getShapeElementAccess().getGroup_3(), "rule__ShapeElement__Group_3__0");
			builder.put(grammarAccess.getPropertyAccess().getGroup(), "rule__Property__Group__0");
			builder.put(grammarAccess.getJsonObjectAccess().getGroup(), "rule__JsonObject__Group__0");
			builder.put(grammarAccess.getJsonObjectAccess().getGroup_1(), "rule__JsonObject__Group_1__0");
			builder.put(grammarAccess.getJsonObjectAccess().getGroup_1_1(), "rule__JsonObject__Group_1_1__0");
			builder.put(grammarAccess.getJsonArrayAccess().getGroup(), "rule__JsonArray__Group__0");
			builder.put(grammarAccess.getJsonArrayAccess().getGroup_1(), "rule__JsonArray__Group_1__0");
			builder.put(grammarAccess.getJsonArrayAccess().getGroup_1_1(), "rule__JsonArray__Group_1_1__0");
			builder.put(grammarAccess.getJsonMemberAccess().getGroup(), "rule__JsonMember__Group__0");
			builder.put(grammarAccess.getLabelElementAccess().getTextAssignment_2_2(), "rule__LabelElement__TextAssignment_2_2");
			builder.put(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_0(), "rule__ElkEdgeSources__SourcesAssignment_1_0");
			builder.put(grammarAccess.getElkEdgeSourcesAccess().getSourcesAssignment_1_1_1(), "rule__ElkEdgeSources__SourcesAssignment_1_1_1");
			builder.put(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_0(), "rule__ElkEdgeTargets__TargetsAssignment_1_0");
			builder.put(grammarAccess.getElkEdgeTargetsAccess().getTargetsAssignment_1_1_1(), "rule__ElkEdgeTargets__TargetsAssignment_1_1_1");
			builder.put(grammarAccess.getElkIdAccess().getIdentifierAssignment_2(), "rule__ElkId__IdentifierAssignment_2");
			builder.put(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_0(), "rule__ElkNodeChildren__ChildrenAssignment_1_0");
			builder.put(grammarAccess.getElkNodeChildrenAccess().getChildrenAssignment_1_1_1(), "rule__ElkNodeChildren__ChildrenAssignment_1_1_1");
			builder.put(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_0(), "rule__ElkNodePorts__PortsAssignment_1_0");
			builder.put(grammarAccess.getElkNodePortsAccess().getPortsAssignment_1_1_1(), "rule__ElkNodePorts__PortsAssignment_1_1_1");
			builder.put(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_0(), "rule__ElkNodeEdges__ContainedEdgesAssignment_1_0");
			builder.put(grammarAccess.getElkNodeEdgesAccess().getContainedEdgesAssignment_1_1_1(), "rule__ElkNodeEdges__ContainedEdgesAssignment_1_1_1");
			builder.put(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_0(), "rule__ElkGraphElementLabels__LabelsAssignment_1_0");
			builder.put(grammarAccess.getElkGraphElementLabelsAccess().getLabelsAssignment_1_1_1(), "rule__ElkGraphElementLabels__LabelsAssignment_1_1_1");
			builder.put(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_0(), "rule__ElkGraphElementProperties__PropertiesAssignment_1_0");
			builder.put(grammarAccess.getElkGraphElementPropertiesAccess().getPropertiesAssignment_1_1_1(), "rule__ElkGraphElementProperties__PropertiesAssignment_1_1_1");
			builder.put(grammarAccess.getShapeElementAccess().getXAssignment_0_2(), "rule__ShapeElement__XAssignment_0_2");
			builder.put(grammarAccess.getShapeElementAccess().getYAssignment_1_2(), "rule__ShapeElement__YAssignment_1_2");
			builder.put(grammarAccess.getShapeElementAccess().getWidthAssignment_2_2(), "rule__ShapeElement__WidthAssignment_2_2");
			builder.put(grammarAccess.getShapeElementAccess().getHeightAssignment_3_2(), "rule__ShapeElement__HeightAssignment_3_2");
			builder.put(grammarAccess.getPropertyAccess().getKeyAssignment_0(), "rule__Property__KeyAssignment_0");
			builder.put(grammarAccess.getPropertyAccess().getValueAssignment_2_0(), "rule__Property__ValueAssignment_2_0");
			builder.put(grammarAccess.getPropertyAccess().getValueAssignment_2_1(), "rule__Property__ValueAssignment_2_1");
			builder.put(grammarAccess.getPropertyAccess().getValueAssignment_2_2(), "rule__Property__ValueAssignment_2_2");
		}
	}
	
	@Inject
	private NameMappings nameMappings;

	@Inject
	private ElkGraphJsonGrammarAccess grammarAccess;

	@Override
	protected InternalElkGraphJsonParser createParser() {
		InternalElkGraphJsonParser result = new InternalElkGraphJsonParser(null);
		result.setGrammarAccess(grammarAccess);
		return result;
	}

	@Override
	protected String getRuleName(AbstractElement element) {
		return nameMappings.getRuleName(element);
	}

	@Override
	protected String[] getInitialHiddenTokens() {
		return new String[] { "RULE_WS", "RULE_ML_COMMENT", "RULE_SL_COMMENT" };
	}

	public ElkGraphJsonGrammarAccess getGrammarAccess() {
		return this.grammarAccess;
	}

	public void setGrammarAccess(ElkGraphJsonGrammarAccess grammarAccess) {
		this.grammarAccess = grammarAccess;
	}
	
	public NameMappings getNameMappings() {
		return nameMappings;
	}
	
	public void setNameMappings(NameMappings nameMappings) {
		this.nameMappings = nameMappings;
	}
}
