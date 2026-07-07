/**
 * ******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
 *  ******************************************************************************
 */
package org.eclipse.elk.core.meta.metaData;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.eclipse.emf.common.util.Enumerator;

/**
 * <!-- begin-user-doc -->
 * A representation of the literals of the enumeration '<em><b>Md Graph Feature</b></em>',
 * and utility methods for working with them.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdGraphFeature()
 * @model
 * @generated
 */
public enum MdGraphFeature implements Enumerator
{
  /**
   * The '<em><b>Self loops</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #SELF_LOOPS_VALUE
   * @generated
   * @ordered
   */
  SELF_LOOPS(0, "self_loops", "self_loops"),

  /**
   * The '<em><b>Inside self loops</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #INSIDE_SELF_LOOPS_VALUE
   * @generated
   * @ordered
   */
  INSIDE_SELF_LOOPS(1, "inside_self_loops", "inside_self_loops"),

  /**
   * The '<em><b>Multi edges</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #MULTI_EDGES_VALUE
   * @generated
   * @ordered
   */
  MULTI_EDGES(2, "multi_edges", "multi_edges"),

  /**
   * The '<em><b>Edge labels</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #EDGE_LABELS_VALUE
   * @generated
   * @ordered
   */
  EDGE_LABELS(3, "edge_labels", "edge_labels"),

  /**
   * The '<em><b>Ports</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #PORTS_VALUE
   * @generated
   * @ordered
   */
  PORTS(4, "ports", "ports"),

  /**
   * The '<em><b>Compound</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #COMPOUND_VALUE
   * @generated
   * @ordered
   */
  COMPOUND(5, "compound", "compound"),

  /**
   * The '<em><b>Clusters</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #CLUSTERS_VALUE
   * @generated
   * @ordered
   */
  CLUSTERS(6, "clusters", "clusters"),

  /**
   * The '<em><b>Disconnected</b></em>' literal object.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #DISCONNECTED_VALUE
   * @generated
   * @ordered
   */
  DISCONNECTED(7, "disconnected", "disconnected");

  /**
   * The '<em><b>Self loops</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Self loops</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #SELF_LOOPS
   * @model name="self_loops"
   * @generated
   * @ordered
   */
  public static final int SELF_LOOPS_VALUE = 0;

  /**
   * The '<em><b>Inside self loops</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Inside self loops</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #INSIDE_SELF_LOOPS
   * @model name="inside_self_loops"
   * @generated
   * @ordered
   */
  public static final int INSIDE_SELF_LOOPS_VALUE = 1;

  /**
   * The '<em><b>Multi edges</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Multi edges</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #MULTI_EDGES
   * @model name="multi_edges"
   * @generated
   * @ordered
   */
  public static final int MULTI_EDGES_VALUE = 2;

  /**
   * The '<em><b>Edge labels</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Edge labels</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #EDGE_LABELS
   * @model name="edge_labels"
   * @generated
   * @ordered
   */
  public static final int EDGE_LABELS_VALUE = 3;

  /**
   * The '<em><b>Ports</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Ports</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #PORTS
   * @model name="ports"
   * @generated
   * @ordered
   */
  public static final int PORTS_VALUE = 4;

  /**
   * The '<em><b>Compound</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Compound</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #COMPOUND
   * @model name="compound"
   * @generated
   * @ordered
   */
  public static final int COMPOUND_VALUE = 5;

  /**
   * The '<em><b>Clusters</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Clusters</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #CLUSTERS
   * @model name="clusters"
   * @generated
   * @ordered
   */
  public static final int CLUSTERS_VALUE = 6;

  /**
   * The '<em><b>Disconnected</b></em>' literal value.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of '<em><b>Disconnected</b></em>' literal object isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @see #DISCONNECTED
   * @model name="disconnected"
   * @generated
   * @ordered
   */
  public static final int DISCONNECTED_VALUE = 7;

  /**
   * An array of all the '<em><b>Md Graph Feature</b></em>' enumerators.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private static final MdGraphFeature[] VALUES_ARRAY =
    new MdGraphFeature[]
    {
      SELF_LOOPS,
      INSIDE_SELF_LOOPS,
      MULTI_EDGES,
      EDGE_LABELS,
      PORTS,
      COMPOUND,
      CLUSTERS,
      DISCONNECTED,
    };

  /**
   * A public read-only list of all the '<em><b>Md Graph Feature</b></em>' enumerators.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public static final List<MdGraphFeature> VALUES = Collections.unmodifiableList(Arrays.asList(VALUES_ARRAY));

  /**
   * Returns the '<em><b>Md Graph Feature</b></em>' literal with the specified literal value.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param literal the literal.
   * @return the matching enumerator or <code>null</code>.
   * @generated
   */
  public static MdGraphFeature get(String literal)
  {
    for (int i = 0; i < VALUES_ARRAY.length; ++i)
    {
      MdGraphFeature result = VALUES_ARRAY[i];
      if (result.toString().equals(literal))
      {
        return result;
      }
    }
    return null;
  }

  /**
   * Returns the '<em><b>Md Graph Feature</b></em>' literal with the specified name.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param name the name.
   * @return the matching enumerator or <code>null</code>.
   * @generated
   */
  public static MdGraphFeature getByName(String name)
  {
    for (int i = 0; i < VALUES_ARRAY.length; ++i)
    {
      MdGraphFeature result = VALUES_ARRAY[i];
      if (result.getName().equals(name))
      {
        return result;
      }
    }
    return null;
  }

  /**
   * Returns the '<em><b>Md Graph Feature</b></em>' literal with the specified integer value.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the integer value.
   * @return the matching enumerator or <code>null</code>.
   * @generated
   */
  public static MdGraphFeature get(int value)
  {
    switch (value)
    {
      case SELF_LOOPS_VALUE: return SELF_LOOPS;
      case INSIDE_SELF_LOOPS_VALUE: return INSIDE_SELF_LOOPS;
      case MULTI_EDGES_VALUE: return MULTI_EDGES;
      case EDGE_LABELS_VALUE: return EDGE_LABELS;
      case PORTS_VALUE: return PORTS;
      case COMPOUND_VALUE: return COMPOUND;
      case CLUSTERS_VALUE: return CLUSTERS;
      case DISCONNECTED_VALUE: return DISCONNECTED;
    }
    return null;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private final int value;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private final String name;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private final String literal;

  /**
   * Only this class can construct instances.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  private MdGraphFeature(int value, String name, String literal)
  {
    this.value = value;
    this.name = name;
    this.literal = literal;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public int getValue()
  {
    return value;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getName()
  {
    return name;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getLiteral()
  {
    return literal;
  }

  /**
   * Returns the literal value of the enumerator, which is its string representation.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public String toString()
  {
    return literal;
  }
  
} //MdGraphFeature
