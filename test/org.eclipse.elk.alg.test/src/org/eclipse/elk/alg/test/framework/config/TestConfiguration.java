/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.config;

import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.StringJoiner;

import org.eclipse.elk.alg.test.framework.annotations.DefaultConfiguration;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.swt.SWT;
import org.junit.runners.model.TestClass;

import com.google.common.base.Strings;

/**
 * Base class of test configurations. This base class knows how to apply defaults and necessary adjustments after layout
 * configurators have been applied. Subclasses must override {@link #applyBaseConfiguration(ElkNode)} and
 * {@link #toString()}.
 */
public class TestConfiguration {
    
    /** Whether to apply default configurations to nodes. */
    private boolean useDefaultNodeConfig = false;
    /** Whether to apply default configurations to ports. */
    private boolean useDefaultPortConfig = false;
    /** Whether to apply default configurations to edges. */
    private boolean useDefaultEdgeConfig = false;
    /** Whether to add random labels to edges. */
    private boolean addRandomEdgeLabels = false;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Creation

    /**
     * Creates a new instance whose defaults are initialized from the given test class's annotations.
     * 
     * @param testClass
     *            test class to initialize this configuration from.
     */
    protected TestConfiguration(final TestClass testClass) {
        DefaultConfiguration dfltConfig = testClass.getAnnotation(DefaultConfiguration.class);
        if (dfltConfig != null) {
            useDefaultNodeConfig = dfltConfig.nodes();
            useDefaultPortConfig = dfltConfig.ports();
            useDefaultEdgeConfig = dfltConfig.edges();
            addRandomEdgeLabels = dfltConfig.edgeLabels();
        }
    }

    /**
     * Returns a test configuration that only applies defaults as configured in the given test class, but leaves all
     * other settings unchanged.
     */
    public static TestConfiguration fromTestClass(final TestClass testClass) {
        return new TestConfiguration(testClass);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TestConfiguration

    /**
     * Applies the strategy to the given graph.
     * 
     * @param test
     *            instance of the test class that runs test methods.
     * @param graph
     *            the graph to apply the strategy to.
     * @throws Throwable
     *             if anything goes wrong.
     */
    public final void apply(final Object test, final ElkNode graph) throws Throwable {
        applyBaseConfiguration(test, graph);
        applyDefaultConfiguration(graph);
        ensureSensibleLabelSize(graph);
    }

    /**
     * Applies basic configuration before defaults are applied.
     * 
     * @param testClass
     *            instance of the test class that runs test methods.
     * @param graph
     *            the graph to configure.
     * @throws Throwable
     *             if anything goes wrong.
     */
    protected void applyBaseConfiguration(final Object testClass, final ElkNode graph) throws Throwable {
        // Default implementation does nothing
    }

    /**
     * Applies defaults to the given graph, if enabled.
     */
    private void applyDefaultConfiguration(final ElkNode graph) {
        if (useDefaultNodeConfig || useDefaultPortConfig || useDefaultEdgeConfig || addRandomEdgeLabels) {
            for (ElkNode elkNode : graph.getChildren()) {
                if (useDefaultNodeConfig) {
                    ElkUtil.configureWithDefaultValues(elkNode);
                }

                if (useDefaultPortConfig || useDefaultEdgeConfig || addRandomEdgeLabels) {
                    for (ElkPort elkPort : elkNode.getPorts()) {
                        if (useDefaultPortConfig) {
                            ElkUtil.configureWithDefaultValues(elkPort);
                        }

                        if (useDefaultEdgeConfig) {
                            for (ElkEdge elkEdge : elkPort.getOutgoingEdges()) {
                                ElkUtil.configureWithDefaultValues(elkEdge);
                            }
                        }
                        
                        if (addRandomEdgeLabels) {
                            for (ElkEdge elkEdge : elkPort.getOutgoingEdges()) {
                                if (elkEdge.getLabels().isEmpty()) {
                                    ElkGraphUtil.createLabel(elkEdge.getIdentifier(), elkEdge);
                                }
                            }
                        }
                    }
                }

                if (useDefaultEdgeConfig) {
                    for (ElkEdge elkEdge : elkNode.getOutgoingEdges()) {
                        ElkUtil.configureWithDefaultValues(elkEdge);
                    }
                }
                
                if (addRandomEdgeLabels) {
                    for (ElkEdge elkEdge : elkNode.getOutgoingEdges()) {
                        if (elkEdge.getLabels().isEmpty()) {
                            ElkGraphUtil.createLabel(elkEdge.getIdentifier(), elkEdge);
                        }
                    }
                }

                applyDefaultConfiguration(elkNode);
            }
        }
    }

    /**
     * Computes a proper size for labels if they have none configured.
     */
    private void ensureSensibleLabelSize(final ElkNode graph) {
        graph.eAllContents().forEachRemaining(obj -> {
            if (obj instanceof ElkLabel) {
                computeLabelSize((ElkLabel) obj);
            }
        });
    }

    /** The size of the font used to calculate the size of a label. */
    private static final int DEFAULT_LABEL_FONT_SIZE = 10;

    /**
     * Computes the size of a label if it doesn't already have one. The size is calculated by choosing a font and
     * calculating the space necessary to render the label's text with that font.
     */
    private void computeLabelSize(final ElkLabel label) {
        // Check if we need to do anything
        if (label.getHeight() == 0 && label.getWidth() == 0 && !Strings.isNullOrEmpty(label.getText())) {
            // There seems to be a problem with SansSerif on Windows...
            String fontName;
            if (System.getProperty("os.name").toLowerCase().indexOf("win") >= 0) {
                fontName = "Arial";
            } else {
                fontName = "SansSerif";
            }

            BufferedImage image = new BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB);

            Graphics2D graphics = image.createGraphics();
            graphics.setFont(new java.awt.Font(fontName, SWT.NORMAL, DEFAULT_LABEL_FONT_SIZE));

            Rectangle2D rect = graphics.getFontMetrics().getStringBounds(label.getText(), graphics);
            label.setHeight(rect.getHeight());
            label.setWidth(rect.getWidth());
        }
    }

    @Override
    public String toString() {
        StringJoiner joiner = new StringJoiner(", ");

        if (useDefaultNodeConfig) {
            joiner.add("nodes");
        }

        if (useDefaultPortConfig) {
            joiner.add("ports");
        }

        if (useDefaultEdgeConfig) {
            joiner.add("edges");
        }
        
        if (addRandomEdgeLabels) {
            joiner.add("edgeLabels");
        }

        return "defaults(" + joiner.toString() + ")";
    }

}
